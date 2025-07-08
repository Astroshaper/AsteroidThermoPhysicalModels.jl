using AsteroidShapeModels
using AsteroidThermoPhysicalModels
using BenchmarkTools
using Downloads
using LinearAlgebra
using Pkg
using Rotations
using SPICE
using StaticArrays
using Statistics

# Import required functions from AsteroidThermoPhysicalModels
using AsteroidThermoPhysicalModels: run_TPM!, update_flux_sun!, update_flux_scat_single!, update_flux_rad_single!, update_temperature!, energy_in, energy_out

# Create benchmark suite
const SUITE = BenchmarkGroup()

# =====================================
# Helper functions for data preparation
# =====================================

"""
Download required SPICE kernels and shape models for benchmarks
"""
function download_benchmark_data()
    # Create directories
    mkpath("benchmark/kernel")
    mkpath("benchmark/shape")
    
    # Ryugu kernels
    ryugu_kernels = [
        "lsk/naif0012.tls",
        "pck/hyb2_ryugu_shape_v20190328.tpc",
        "fk/hyb2_ryugu_v01.tf",
        "spk/2162173_Ryugu.bsp",
    ]
    
    for path_kernel in ryugu_kernels
        url_kernel = "https://data.darts.isas.jaxa.jp/pub/hayabusa2/old/2020/spice_bundle/spice_kernels/$(path_kernel)"
        filepath = joinpath("benchmark/kernel", path_kernel)
        mkpath(dirname(filepath))
        if !isfile(filepath)
            @info "Downloading $(path_kernel)..."
            Downloads.download(url_kernel, filepath)
        end
    end
    
    # Ryugu shape model
    ryugu_shapes = ["SHAPE_SFM_49k_v20180804.obj"]
    for path_shape in ryugu_shapes
        url_shape = "https://data.darts.isas.jaxa.jp/pub/hayabusa2/paper/Watanabe_2019/$(path_shape)"
        filepath = joinpath("benchmark/shape", path_shape)
        if !isfile(filepath)
            @info "Downloading $(path_shape)..."
            Downloads.download(url_shape, filepath)
        end
    end
    
    # Didymos kernels
    didymos_kernels = [
        "fk/hera_v10.tf",
        "lsk/naif0012.tls",
        "pck/hera_didymos_v06.tpc",
        "spk/de432s.bsp",
        "spk/didymos_hor_000101_500101_v01.bsp",
        "spk/didymos_gmv_260901_311001_v01.bsp",
    ]
    
    for path_kernel in didymos_kernels
        url_kernel = "https://s2e2.cosmos.esa.int/bitbucket/projects/SPICE_KERNELS/repos/hera/raw/kernels/$(path_kernel)?at=refs%2Ftags%2Fv161_20230929_001"
        filepath = joinpath("benchmark/kernel", path_kernel)
        mkpath(dirname(filepath))
        if !isfile(filepath)
            @info "Downloading $(path_kernel)..."
            Downloads.download(url_kernel, filepath)
        end
    end
    
    # Didymos shape models
    didymos_shapes = [
        "g_50677mm_rad_obj_didy_0000n00000_v001.obj",
        "g_08438mm_lgt_obj_dimo_0000n00000_v002.obj",
    ]
    
    for path_shape in didymos_shapes
        url_kernel = "https://s2e2.cosmos.esa.int/bitbucket/projects/SPICE_KERNELS/repos/hera/raw/kernels/dsk/$(path_shape)?at=refs%2Ftags%2Fv161_20230929_001"
        filepath = joinpath("benchmark/shape", path_shape)
        if !isfile(filepath)
            @info "Downloading $(path_shape)..."
            Downloads.download(url_kernel, filepath)
        end
    end
end

"""
Setup Ryugu TPM for benchmarking
"""
function setup_ryugu_tpm(n_cycles::Int=20)
    # Load SPICE kernels
    kernel_paths = [
        "lsk/naif0012.tls",
        "pck/hyb2_ryugu_shape_v20190328.tpc",
        "fk/hyb2_ryugu_v01.tf",
        "spk/2162173_Ryugu.bsp",
    ]
    
    for path_kernel in kernel_paths
        filepath = joinpath("benchmark/kernel", path_kernel)
        SPICE.furnsh(filepath)
    end
    
    # Ephemerides
    P = SPICE.convrt(7.63262, "hours", "seconds")  # Rotation period of Ryugu
    n_step_in_cycle = 72  # Number of time steps in one rotation period
    
    et_begin = SPICE.utc2et("2018-07-01T00:00:00")
    et_end = et_begin + P * n_cycles
    et_range = range(et_begin, et_end; length=n_step_in_cycle*n_cycles+1)
    
    ephem = (
        time = collect(et_range),
        sun = [SVector{3}(SPICE.spkpos("SUN", et, "RYUGU_FIXED", "None", "RYUGU")[1]) * 1000 for et in et_range],
    )
    
    SPICE.kclear()
    
    # Load shape model
    path_obj = joinpath("benchmark/shape", "SHAPE_SFM_49k_v20180804.obj")
    shape = load_shape_obj(path_obj; scale=1000, with_face_visibility=true)
    
    # Thermal properties
    k = 0.1      # Thermal conductivity [W/m/K]
    ρ = 1270.0   # Density [kg/m³]
    Cₚ = 600.0   # Heat capacity [J/kg/K]
    
    l = AsteroidThermoPhysicalModels.thermal_skin_depth(P, k, ρ, Cₚ)
    Γ = AsteroidThermoPhysicalModels.thermal_inertia(k, ρ, Cₚ)
    
    R_vis = 0.04  # Reflectance in visible light
    R_ir = 0.0    # Reflectance in thermal infrared
    ε = 1.0       # Emissivity
    
    z_max = 0.6   # Depth of the lower boundary [m]
    n_depth = 41  # Number of depth steps
    Δz = z_max / (n_depth - 1)
    
    thermo_params = AsteroidThermoPhysicalModels.ThermoParams(k, ρ, Cₚ, R_vis, R_ir, ε, z_max, Δz, n_depth)
    
    # Create TPM with all features enabled for realistic benchmark
    stpm = AsteroidThermoPhysicalModels.SingleAsteroidTPM(shape, thermo_params;
        SELF_SHADOWING = true,   # Enable shadow calculations
        SELF_HEATING = true,     # Enable self-heating
        SOLVER = AsteroidThermoPhysicalModels.CrankNicolsonSolver(thermo_params),
        BC_UPPER = AsteroidThermoPhysicalModels.RadiationBoundaryCondition(),
        BC_LOWER = AsteroidThermoPhysicalModels.InsulationBoundaryCondition(),
    )
    
    AsteroidThermoPhysicalModels.init_temperature!(stpm, 200)
    
    return stpm, ephem
end

"""
Setup Didymos-Dimorphos binary TPM for benchmarking
"""
function setup_didymos_tpm(n_cycles::Int=20)
    # Load SPICE kernels
    kernel_paths = [
        "fk/hera_v10.tf",
        "lsk/naif0012.tls",
        "pck/hera_didymos_v06.tpc",
        "spk/de432s.bsp",
        "spk/didymos_hor_000101_500101_v01.bsp",
        "spk/didymos_gmv_260901_311001_v01.bsp",
    ]
    
    for path_kernel in kernel_paths
        filepath = joinpath("benchmark/kernel", path_kernel)
        SPICE.furnsh(filepath)
    end
    
    # Ephemerides
    P₁ = SPICE.convrt(2.2593, "hours", "seconds")  # Rotation period of Didymos
    P₂ = SPICE.convrt(11.93, "hours", "seconds")   # Rotation period of Dimorphos
    
    n_step_in_cycle = 72
    
    et_begin = SPICE.utc2et("2027-02-18T00:00:00")
    et_end = et_begin + P₂ * n_cycles  # Use secondary's period
    et_range = range(et_begin, et_end; length=n_step_in_cycle*n_cycles+1)
    
    ephem = (
        time = collect(et_range),
        sun1 = [SVector{3}(SPICE.spkpos("SUN", et, "DIDYMOS_FIXED", "None", "DIDYMOS")[1]) * 1000 for et in et_range],
        sun2 = [SVector{3}(SPICE.spkpos("SUN", et, "DIMORPHOS_FIXED", "None", "DIMORPHOS")[1]) * 1000 for et in et_range],
        sec = [SVector{3}(SPICE.spkpos("DIMORPHOS", et, "DIDYMOS_FIXED", "None", "DIDYMOS")[1]) * 1000 for et in et_range],
        P2S = [RotMatrix{3}(SPICE.pxform("DIDYMOS_FIXED", "DIMORPHOS_FIXED", et)) for et in et_range],
        S2P = [RotMatrix{3}(SPICE.pxform("DIMORPHOS_FIXED", "DIDYMOS_FIXED", et)) for et in et_range],
    )
    
    SPICE.kclear()
    
    # Load shape models
    path_shape1_obj = joinpath("benchmark/shape", "g_50677mm_rad_obj_didy_0000n00000_v001.obj")
    path_shape2_obj = joinpath("benchmark/shape", "g_08438mm_lgt_obj_dimo_0000n00000_v002.obj")
    
    shape1 = load_shape_obj(path_shape1_obj; scale=1000, with_face_visibility=true)
    shape2 = load_shape_obj(path_shape2_obj; scale=1000, with_face_visibility=true)
    
    # Thermal properties
    k = 0.125    # Thermal conductivity [W/m/K]
    ρ = 2170.0   # Density [kg/m³]
    Cₚ = 600.0   # Heat capacity [J/kg/K]
    
    l₁ = AsteroidThermoPhysicalModels.thermal_skin_depth(P₁, k, ρ, Cₚ)
    l₂ = AsteroidThermoPhysicalModels.thermal_skin_depth(P₂, k, ρ, Cₚ)
    Γ = AsteroidThermoPhysicalModels.thermal_inertia(k, ρ, Cₚ)
    
    R_vis = 0.059  # Reflectance in visible light
    R_ir = 0.0     # Reflectance in thermal infrared
    ε = 0.9        # Emissivity
    
    z_max = 0.6
    n_depth = 41
    Δz = z_max / (n_depth - 1)
    
    thermo_params1 = AsteroidThermoPhysicalModels.ThermoParams(k, ρ, Cₚ, R_vis, R_ir, ε, z_max, Δz, n_depth)
    thermo_params2 = AsteroidThermoPhysicalModels.ThermoParams(k, ρ, Cₚ, R_vis, R_ir, ε, z_max, Δz, n_depth)
    
    # Create TPMs with all features enabled
    stpm1 = AsteroidThermoPhysicalModels.SingleAsteroidTPM(shape1, thermo_params1;
        SELF_SHADOWING = true,
        SELF_HEATING = true,
        SOLVER = AsteroidThermoPhysicalModels.CrankNicolsonSolver(thermo_params1),
        BC_UPPER = AsteroidThermoPhysicalModels.RadiationBoundaryCondition(),
        BC_LOWER = AsteroidThermoPhysicalModels.InsulationBoundaryCondition(),
    )
    
    stpm2 = AsteroidThermoPhysicalModels.SingleAsteroidTPM(shape2, thermo_params2;
        SELF_SHADOWING = true,
        SELF_HEATING = true,
        SOLVER = AsteroidThermoPhysicalModels.CrankNicolsonSolver(thermo_params2),
        BC_UPPER = AsteroidThermoPhysicalModels.RadiationBoundaryCondition(),
        BC_LOWER = AsteroidThermoPhysicalModels.InsulationBoundaryCondition(),
    )
    
    btpm = AsteroidThermoPhysicalModels.BinaryAsteroidTPM(stpm1, stpm2;
        MUTUAL_SHADOWING = true,  # Enable mutual eclipses
        MUTUAL_HEATING = true,    # Enable mutual heating
    )
    
    AsteroidThermoPhysicalModels.init_temperature!(btpm, 200.0)
    
    return btpm, ephem
end

# =====================================
# Benchmark: Ryugu simulation
# =====================================

SUITE["ryugu"] = BenchmarkGroup()

# Full 20-rotation simulation
SUITE["ryugu"]["full_simulation_20_rotations"] = @benchmarkable begin
    result = run_TPM!(stpm, ephem, times_to_save, face_ID; show_progress=false)
end setup = begin
    @info "Setting up Ryugu benchmark (20 rotations)..."
    stpm, ephem = setup_ryugu_tpm(20)
    times_to_save = Float64[]  # Don't save intermediate results for benchmark
    face_ID = Int[]
end teardown = begin
    # Cleanup if needed
end

# Single rotation for comparison
SUITE["ryugu"]["single_rotation"] = @benchmarkable begin
    result = run_TPM!(stpm, ephem, times_to_save, face_ID; show_progress=false)
end setup = begin
    @info "Setting up Ryugu benchmark (1 rotation)..."
    stpm, ephem = setup_ryugu_tpm(1)
    times_to_save = Float64[]
    face_ID = Int[]
end

# =====================================
# Benchmark: Didymos-Dimorphos binary
# =====================================

SUITE["didymos"] = BenchmarkGroup()

# Full 20-rotation simulation
SUITE["didymos"]["full_simulation_20_rotations"] = @benchmarkable begin
    result = run_TPM!(btpm, ephem, times_to_save, face_ID_pri, face_ID_sec; show_progress=false)
end setup = begin
    @info "Setting up Didymos-Dimorphos benchmark (20 rotations)..."
    btpm, ephem = setup_didymos_tpm(20)
    times_to_save = Float64[]
    face_ID_pri = Int[]
    face_ID_sec = Int[]
end teardown = begin
    # Cleanup if needed
end

# Single rotation for comparison
SUITE["didymos"]["single_rotation"] = @benchmarkable begin
    result = run_TPM!(btpm, ephem, times_to_save, face_ID_pri, face_ID_sec; show_progress=false)
end setup = begin
    @info "Setting up Didymos-Dimorphos benchmark (1 rotation)..."
    btpm, ephem = setup_didymos_tpm(1)
    times_to_save = Float64[]
    face_ID_pri = Int[]
    face_ID_sec = Int[]
end

# =====================================
# Benchmark: Component analysis
# =====================================

SUITE["components"] = BenchmarkGroup()

# Shadow calculation overhead
SUITE["components"]["ryugu_shadow_overhead"] = @benchmarkable begin
    for i in 1:length(ephem.time)
        update_flux_sun!(stpm, ephem.sun[i])
    end
end setup = begin
    stpm, ephem = setup_ryugu_tpm(1)
end

# Self-heating calculation overhead
SUITE["components"]["ryugu_self_heating"] = @benchmarkable begin
    for i in 1:length(ephem.time)
        update_flux_scat_single!(stpm)
        update_flux_rad_single!(stpm)
    end
end setup = begin
    stpm, ephem = setup_ryugu_tpm(1)
    # Initialize with some flux values
    stpm.flux_sun .= 1000.0
    AsteroidThermoPhysicalModels.init_temperature!(stpm, 250.0)
end

# Temperature update performance
SUITE["components"]["ryugu_temperature_update"] = @benchmarkable begin
    for i in 1:length(ephem.time)-1
        Δt = ephem.time[i+1] - ephem.time[i]
        update_temperature!(stpm, Δt)
    end
end setup = begin
    stpm, ephem = setup_ryugu_tpm(1)
    stpm.flux_sun .= 1000.0
end

# =====================================
# Memory benchmarks
# =====================================

SUITE["memory"] = BenchmarkGroup()

# Memory allocation in full simulation
SUITE["memory"]["ryugu_full_simulation"] = @benchmarkable begin
    result = run_TPM!(stpm, ephem, times_to_save, face_ID; show_progress=false)
end setup = begin
    stpm, ephem = setup_ryugu_tpm(1)
    times_to_save = ephem.time[end:end]
    face_ID = [1, 100, 1000]
end

# =====================================
# Utility functions
# =====================================

"""
Run all benchmarks and save results
"""
function run_benchmarks(; save_results=true)
    # Download data if needed
    download_benchmark_data()
    
    # Run benchmarks
    results = run(SUITE, verbose=true)
    
    if save_results
        # Get package version
        pkg_version = Pkg.project().version
        
        # Save results with timestamp and version
        timestamp = Dates.format(now(), "yyyy-mm-dd_HHMMSS")
        filename = "benchmark_results_v$(pkg_version)_$(timestamp).json"
        BenchmarkTools.save(filename, results)
        @info "Benchmark results saved to $(filename)"
    end
    
    return results
end

"""
Compare benchmark results between two runs
"""
function compare_benchmarks(file1::String, file2::String)
    results1 = BenchmarkTools.load(file1)
    results2 = BenchmarkTools.load(file2)
    
    comparison = judge(results2, results1)
    
    return comparison
end

"""
Print benchmark summary
"""
function print_summary(results)
    println("\n" * "="^60)
    println("BENCHMARK SUMMARY")
    println("="^60)
    
    # Ryugu results
    if haskey(results, "ryugu")
        println("\nRyugu Simulation (49,152 faces):")
        if haskey(results["ryugu"], "single_rotation")
            r1 = results["ryugu"]["single_rotation"]
            t1 = median(r1).time / 1e9  # Convert to seconds
            println("  1 rotation (72 steps):")
            println("    Time        : $(BenchmarkTools.prettytime(median(r1).time))")
            println("    Memory      : $(BenchmarkTools.prettymemory(median(r1).memory))")
            println("    Allocations : $(median(r1).allocs)")
        end
        if haskey(results["ryugu"], "full_simulation_20_rotations")
            r20 = results["ryugu"]["full_simulation_20_rotations"]
            t20 = median(r20).time / 1e9  # Convert to seconds
            println("  - 20 rotations (1,440 steps):")
            println("      - Time        : $(BenchmarkTools.prettytime(median(r20).time))")
            println("      - Memory      : $(BenchmarkTools.prettymemory(median(r20).memory))")
            println("      - Allocations : $(median(r20).allocs)")
        end
    end
    
    # Didymos results
    if haskey(results, "didymos")
        println("\nDidymos-Dimorphos Binary Simulation (1,996 + 3,072 faces):")
        if haskey(results["didymos"], "single_rotation")
            r1 = results["didymos"]["single_rotation"]
            t1 = median(r1).time / 1e9  # Convert to seconds
            println("  - 1 rotation (72 steps):")
            println("      - Time        : $(BenchmarkTools.prettytime(median(r1).time))")
            println("      - Memory      : $(BenchmarkTools.prettymemory(median(r1).memory))")
            println("      - Allocations : $(median(r1).allocs)")
        end
        if haskey(results["didymos"], "full_simulation_20_rotations")
            r20 = results["didymos"]["full_simulation_20_rotations"]
            t20 = median(r20).time / 1e9  # Convert to seconds
            println("  - 20 rotations (1,440 steps):")
            println("      - Time        : $(BenchmarkTools.prettytime(median(r20).time))")
            println("      - Memory      : $(BenchmarkTools.prettymemory(median(r20).memory))")
            println("      - Allocations : $(median(r20).allocs)")
        end
    end
    
    println("\n" * "="^60)
end
