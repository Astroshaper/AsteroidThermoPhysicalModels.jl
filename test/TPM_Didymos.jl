# See https://github.com/Astroshaper/Astroshaper-examples/tree/main/TPM_Didymos for more information.
@testset "TPM_Didymos" begin
    msg = """\n
    ┏━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━┓
    |                   Test: TPM_Didymos                    |
    ┗━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━┛
    """
    println(msg)

    ##= SPICE kernels =##
    paths_kernel = [
        "fk/hera_v10.tf",
        "lsk/naif0012.tls",
        "pck/hera_didymos_v06.tpc",
        "spk/de432s.bsp",
        "spk/didymos_hor_000101_500101_v01.bsp",
        "spk/didymos_gmv_260901_311001_v01.bsp",
    ]

    ##= Shape models =##
    paths_shape = [
        "g_50677mm_rad_obj_didy_0000n00000_v001.obj",
        "g_08438mm_lgt_obj_dimo_0000n00000_v002.obj",
    ]

    ##= Download SPICE kernels =##
    for path_kernel in paths_kernel
        url_kernel = "https://s2e2.cosmos.esa.int/bitbucket/projects/SPICE_KERNELS/repos/hera/raw/kernels/$(path_kernel)"
        filepath = joinpath("kernel", path_kernel)
        mkpath(dirname(filepath))
        isfile(filepath) || Downloads.download(url_kernel, filepath)
    end

    ##= Download shape models =##
    for path_shape in paths_shape
        url_kernel = "https://s2e2.cosmos.esa.int/bitbucket/projects/SPICE_KERNELS/repos/hera/raw/kernels/dsk/$(path_shape)"
        filepath = joinpath("shape", path_shape)
        mkpath(dirname(filepath))
        isfile(filepath) || Downloads.download(url_kernel, filepath)
    end

    ##= Load the SPICE kernels =##
    for path_kernel in paths_kernel
        filepath = joinpath("kernel", path_kernel)
        SPICE.furnsh(filepath)
    end

    ##= Ephemerides =##
    P₁ = SPICE.convrt(2.2593, "hours", "seconds")  # Rotation period of Didymos
    P₂ = SPICE.convrt(11.93 , "hours", "seconds")  # Rotation period of Dimorphos

    et_begin = SPICE.utc2et("2027-02-18T00:00:00")  # Start time of TPM
    et_end   = et_begin + 2P₂                       # End time of TPM

    nsteps_in_rotation = 72  # Number of steps in one rotation period
    step = P₂ / nsteps_in_rotation  # Time step of TPM
    et_range = et_begin : step : et_end

    """
    - `time` : Ephemeris times
    - `sun1` : Sun's position in the primary's frame (DIDYMOS_FIXED)
    - `sun2` : Sun's position in the secondary's frame (DIMORPHOS_FIXED)
    - `sec`  : Secondary's position in the primary's frame (DIDYMOS_FIXED)
    - `P2S`  : Rotation matrix from primary to secondary frames
    - `S2P`  : Rotation matrix from secondary to primary frames
    """
    ephem = (
        time = collect(et_range),
        sun1 = [SVector{3}(SPICE.spkpos("SUN"      , et, "DIDYMOS_FIXED"  , "None", "DIDYMOS"  )[1]) * 1000 for et in et_range],
        sun2 = [SVector{3}(SPICE.spkpos("SUN"      , et, "DIMORPHOS_FIXED", "None", "DIMORPHOS")[1]) * 1000 for et in et_range],
        sec  = [SVector{3}(SPICE.spkpos("DIMORPHOS", et, "DIDYMOS_FIXED"  , "None", "DIDYMOS"  )[1]) * 1000 for et in et_range],
        P2S  = [RotMatrix{3}(SPICE.pxform("DIDYMOS_FIXED"  , "DIMORPHOS_FIXED", et)) for et in et_range],
        S2P  = [RotMatrix{3}(SPICE.pxform("DIMORPHOS_FIXED", "DIDYMOS_FIXED"  , et)) for et in et_range],
    )

    SPICE.kclear()

    ##= Load the shape models =##
    path_shape1_obj = joinpath("shape", "g_50677mm_rad_obj_didy_0000n00000_v001.obj")
    path_shape2_obj = joinpath("shape", "g_08438mm_lgt_obj_dimo_0000n00000_v002.obj")
    
    shape1 = AsteroidThermoPhysicalModels.load_shape_obj(path_shape1_obj; scale=1000, find_visible_facets=true)
    shape2 = AsteroidThermoPhysicalModels.load_shape_obj(path_shape2_obj; scale=1000, find_visible_facets=true)
    
    ##= Thermal properties =##
    k  = 0.125
    ρ  = 2170.
    Cₚ = 600.

    l₁ = AsteroidThermoPhysicalModels.thermal_skin_depth(P₁, k, ρ, Cₚ)  # Thermal skin depth for Didymos
    l₂ = AsteroidThermoPhysicalModels.thermal_skin_depth(P₂, k, ρ, Cₚ)  # Thermal skin depth for Dimorphos
    Γ = AsteroidThermoPhysicalModels.thermal_inertia(k, ρ, Cₚ)

    thermo_params1 = AsteroidThermoPhysicalModels.thermoparams(  # [Michel+2016; Naidu+2020]
        P       = P₁,
        l       = l₁,
        Γ       = Γ,
        A_B     = 0.059,  # Bolometric Bond albedo
        A_TH    = 0.0,
        ε       = 0.9,
        z_max   = 0.6,
        Nz      = 41,
    )

    thermo_params2 = AsteroidThermoPhysicalModels.thermoparams(  # [Michel+2016; Naidu+2020]
        P       = P₂,
        l       = l₂,
        Γ       = Γ,
        A_B     = 0.059,  # Bolometric Bond albedo
        A_TH    = 0.0,
        ε       = 0.9,
        z_max   = 0.6,
        Nz      = 41,
    )

    println("Thermophysical parameters for Didymos")
    println(thermo_params1)
    println("Thermophysical parameters for Dimorphos")
    println(thermo_params2)

    ##= Setting of TPM =##
    stpm1 = AsteroidThermoPhysicalModels.SingleTPM(shape1, thermo_params1;
        SELF_SHADOWING = true,
        SELF_HEATING   = true,
        SOLVER         = AsteroidThermoPhysicalModels.ForwardEulerSolver(thermo_params1),
        BC_UPPER       = AsteroidThermoPhysicalModels.RadiationBoundaryCondition(),
        BC_LOWER       = AsteroidThermoPhysicalModels.InsulationBoundaryCondition(),
    )

    stpm2 = AsteroidThermoPhysicalModels.SingleTPM(shape2, thermo_params2;
        SELF_SHADOWING = true,
        SELF_HEATING   = true,
        SOLVER         = AsteroidThermoPhysicalModels.ForwardEulerSolver(thermo_params2),
        BC_UPPER       = AsteroidThermoPhysicalModels.RadiationBoundaryCondition(),
        BC_LOWER       = AsteroidThermoPhysicalModels.InsulationBoundaryCondition(),
    )

    btpm  = AsteroidThermoPhysicalModels.BinaryTPM(stpm1, stpm2; MUTUAL_SHADOWING=true, MUTUAL_HEATING=true)
    AsteroidThermoPhysicalModels.init_temperature!(btpm, 200.)
    
    ##= Run TPM =##
    times_to_save = ephem.time[end-nsteps_in_rotation:end]  # Save temperature during the final rotation
    face_ID_pri = [1, 2, 3, 4, 10]  # Face indices to save subsurface temperature of the primary
    face_ID_sec = [1, 2, 3, 4, 20]  # Face indices to save subsurface temperature of the secondary

    result = AsteroidThermoPhysicalModels.run_TPM!(btpm, ephem, times_to_save, face_ID_pri, face_ID_sec)

    ##= Save TPM result =##
    savedir = "TPM_Didymos"
    mkpath(savedir)
    AsteroidThermoPhysicalModels.export_TPM_results(savedir, result, btpm, ephem, times_to_save)
end
