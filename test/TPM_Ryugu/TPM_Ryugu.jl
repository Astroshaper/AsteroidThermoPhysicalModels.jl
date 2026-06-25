#=
TPM_Ryugu.jl

End-to-end thermophysical model test for asteroid Ryugu.
This test simulates:
- Temperature evolution using actual Ryugu shape model
- Self-shadowing and self-heating effects
- Validation of TPM implementation for real asteroids
See https://github.com/Astroshaper/Astroshaper-examples/tree/main/TPM_Ryugu for more information.
=#

@testset "TPM_Ryugu" begin
    DIR_OUTPUT = joinpath(@__DIR__, "output")
    rm(DIR_OUTPUT; recursive=true, force=true)
    mkpath(DIR_OUTPUT)

    msg = """\n
    ┏━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━┓
    |                    Test: TPM_Ryugu                     |
    ┗━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━┛
    """
    println(msg)

    ## --- Download Files ---
    paths_kernel = [
        "lsk/naif0012.tls",
        "pck/hyb2_ryugu_shape_v20190328.tpc",
        "fk/hyb2_ryugu_v01.tf",
        "spk/2162173_Ryugu.bsp",
    ]
    paths_shape = [
        "SHAPE_SFM_49k_v20180804.obj",
    ]

    for path_kernel in paths_kernel
        url_kernel = "https://data.darts.isas.jaxa.jp/pub/hayabusa2/old/2020/spice_bundle/spice_kernels/$(path_kernel)"
        filepath = joinpath("kernel", path_kernel)
        mkpath(dirname(filepath))
        isfile(filepath) || Downloads.download(url_kernel, filepath)
    end

    for path_shape in paths_shape
        url_shape = "https://data.darts.isas.jaxa.jp/pub/hayabusa2/paper/Watanabe_2019/$(path_shape)"
        filepath = joinpath("shape", path_shape)
        mkpath(dirname(filepath))
        isfile(filepath) || Downloads.download(url_shape, filepath)
    end

    ## --- Load data with SPICE ---
    for path_kernel in paths_kernel
        filepath = joinpath("kernel", path_kernel)
        SPICE.furnsh(filepath)
    end

    ## --- Ephemerides ---
    P = SPICE.convrt(7.63262, "hours", "seconds")  # Rotation period of Ryugu

    n_cycle = 2  # Number of cycles to perform TPM
    n_step_in_cycle = 120  # Number of time steps in one rotation period

    et_begin = SPICE.utc2et("2018-07-01T00:00:00")  # Start time of TPM
    et_end   = et_begin + P * n_cycle  # End time of TPM
    et_range = range(et_begin, et_end; length=n_step_in_cycle*n_cycle+1)

    times = collect(et_range)
    r_sun = [SPICE.spkpos("SUN", et, "RYUGU_FIXED", "None", "RYUGU")[1] for et in et_range]  # Sun's position in RYUGU_FIXED [km]
    r_sun .*= 1000  # Convert [km] to [m]


    ephem = SingleAsteroidEphemerides(times, r_sun)

    SPICE.kclear()

    ## --- Load obj file ---
    path_obj = joinpath("shape", "ryugu_test.obj")  # Small model for test
    # path_obj = joinpath("shape", "SHAPE_SFM_49k_v20180804.obj")
    
    shape = load_shape_obj(path_obj; scale=1000, with_face_visibility=true, with_bvh=true)
    n_face = length(shape.faces)  # Number of faces

    ## --- Thermal properties ---
    k  = 0.1     # Thermal conductivity [W/m/K]
    ρ  = 1270.0  # Density [kg/m³]
    Cₚ = 600.0   # Heat capacity [J/kg/K]

    l = thermal_skin_depth(P, k, ρ, Cₚ)  # Thermal skin depth [m]
    Γ = thermal_inertia(k, ρ, Cₚ)        # Thermal inertia [tiu]

    R_vis = 0.04  # Reflectance in visible light [-]
    R_ir  = 0.0   # Reflectance in thermal infrared [-]
    ε     = 1.0   # Emissivity [-]

    z_max = 0.6   # Depth of the lower boundary of a heat conduction equation [m]
    n_depth = 61  # Number of depth steps
    Δz = z_max / (n_depth - 1)  # Depth step width [m]

    thermo_params = ThermoParams(k, ρ, Cₚ, R_vis, R_ir, ε, z_max, Δz, n_depth)

    ## --- Setting of TPM ---
    problem = SingleAsteroidThermoPhysicalProblem(shape, thermo_params;
        with_self_shadowing      = true,
        with_self_heating        = true,
        upper_boundary_condition = RadiationBoundaryCondition(),
        lower_boundary_condition = InsulationBoundaryCondition(),
    )

    ## --- Run TPM ---
    output_times        = ephem.times[end-n_step_in_cycle:end]  # Save temperature during the final rotation
    subsurface_face_ids = [1, 2, 3, 4, 10]  # Face indices to save subsurface temperature
    output = SingleAsteroidOutputSpec(output_times, subsurface_face_ids)

    solution = solve(problem, ExplicitEuler();
        ephem               = ephem,
        output              = output,
        initial_temperature = 200.0,
    )

    ## --- Save TPM result ---
    @testset "Save TPM result" begin
        export_solution(DIR_OUTPUT, solution)

        @test isfile(joinpath(DIR_OUTPUT, "diagnostics.csv"))
        @test isfile(joinpath(DIR_OUTPUT, "surface_temperature.csv"))
        @test isfile(joinpath(DIR_OUTPUT, "subsurface_temperature.csv"))
        @test !isfile(joinpath(DIR_OUTPUT, "thermal_face_forces.csv"))  # save_face_forces=false by default
    end
end
