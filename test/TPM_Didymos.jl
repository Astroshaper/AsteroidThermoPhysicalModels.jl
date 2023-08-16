# See https://github.com/Astroshaper/Astroshaper-examples/tree/main/TPM_Didymos for more information.
@testset "TPM_Didymos" begin
    msg = """\n
    ┏━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━┓
    |                   Test: TPM_Didymos                    |
    ┗━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━┛
    """
    println(msg)

    ##= Download Files =##
    paths_kernel = [
        "fk/hera_v09.tf",
        "lsk/naif0012.tls",
        "pck/hera_didymos_v06.tpc",
        "spk/de432s.bsp",
        "spk/didymos_hor_000101_500101_v01.bsp",
        "spk/didymos_gmv_260901_311001_v01.bsp",
    ]
    paths_shape = [
        "g_50677mm_rad_obj_didy_0000n00000_v001.obj",
        "g_08438mm_lgt_obj_dimo_0000n00000_v002.obj",
    ]

    for path_kernel in paths_kernel
        url_kernel = "https://s2e2.cosmos.esa.int/bitbucket/projects/SPICE_KERNELS/repos/hera/raw/kernels/$(path_kernel)"
        filepath = joinpath("kernel", path_kernel)
        mkpath(dirname(filepath))
        isfile(filepath) || Downloads.download(url_kernel, filepath)
    end
    for path_shape in paths_shape
        url_kernel = "https://s2e2.cosmos.esa.int/bitbucket/projects/SPICE_KERNELS/repos/hera/raw/kernels/dsk/$(path_shape)"
        filepath = joinpath("shape", path_shape)
        mkpath(dirname(filepath))
        isfile(filepath) || Downloads.download(url_kernel, filepath)
    end

    ##= Load data with SPICE =##
    for path_kernel in paths_kernel
        filepath = joinpath("kernel", path_kernel)
        SPICE.furnsh(filepath)
    end

    ##= Ephemerides =##
    P        = SPICE.convrt(AsteroidThermoPhysicalModels.DIDYMOS[:P], "hours", "seconds")  # Rotation period of Didymos
    et_begin = SPICE.utc2et("2027-02-18T00:00:00")                                         # Start time of TPM
    et_end   = et_begin + 10P                                                              # End time of TPM
    step     = P / 72                                                                      # Time step of TPM
    et_range = et_begin : step : et_end
    @show length(et_range)

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

    ##= Load obj file =##
    path_shape1_obj = joinpath("shape", "g_50677mm_rad_obj_didy_0000n00000_v001.obj")
    path_shape2_obj = joinpath("shape", "g_08438mm_lgt_obj_dimo_0000n00000_v002.obj")
    path_shape1_jld = joinpath("shape", "g_50677mm_rad_obj_didy_0000n00000_v001.jld2")
    path_shape2_jld = joinpath("shape", "g_08438mm_lgt_obj_dimo_0000n00000_v002.jld2")

    if isfile(path_shape1_jld) && ENABLE_JLD
        shape1 = AsteroidThermoPhysicalModels.load_shape_jld(path_shape1_jld)
    else
        shape1 = AsteroidThermoPhysicalModels.load_shape_obj(path_shape1_obj; scale=1000, find_visible_facets=true)
        AsteroidThermoPhysicalModels.save_shape_jld(path_shape1_jld, shape1)
    end
    if isfile(path_shape2_jld) && ENABLE_JLD
        shape2 = AsteroidThermoPhysicalModels.load_shape_jld(path_shape2_jld)
    else
        shape2 = AsteroidThermoPhysicalModels.load_shape_obj(path_shape2_obj; scale=1000, find_visible_facets=true)
        AsteroidThermoPhysicalModels.save_shape_jld(path_shape2_jld, shape2)
    end
    
    ##= Thermal properties =##
    P  = SPICE.convrt(AsteroidThermoPhysicalModels.DIDYMOS[:P], "hours", "seconds")
    k  = 0.125
    ρ  = 2170.
    Cₚ = 600.

    l = AsteroidThermoPhysicalModels.thermal_skin_depth(P, k, ρ, Cₚ)
    Γ = AsteroidThermoPhysicalModels.thermal_inertia(k, ρ, Cₚ)

    thermo_params = AsteroidThermoPhysicalModels.thermoparams(  # [Michel+2016; Naidu+2020]
        P       = P,
        l       = l,
        Γ       = Γ,
        A_B     = 0.059,  # Bolometric Bond albedo
        A_TH    = 0.0,
        ε       = 0.9,
        t_begin = et_range[begin],
        t_end   = et_range[end],
        Nt      = length(et_range),
        z_max   = 0.6,
        Nz      = 41,
    )

    println(thermo_params)

    ##= Setting of TPM =##
    tpm1 = AsteroidThermoPhysicalModels.SingleTPM(shape1, thermo_params, true, true)
    tpm2 = AsteroidThermoPhysicalModels.SingleTPM(shape2, thermo_params, true, true)
    btpm = AsteroidThermoPhysicalModels.BinaryTPM(tpm1, tpm2, true, true)

    AsteroidThermoPhysicalModels.init_temperature!(btpm, 200.)

    ##= Run TPM and save the result =##
    savepath = "TPM_Didymos.jld2"
    AsteroidThermoPhysicalModels.run_TPM!(btpm, ephem, savepath)
end
