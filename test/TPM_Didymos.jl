# See https://github.com/Astroshaper/Astroshaper-examples/tree/main/TPM_Didymos for more information.
@testset "TPM_Didymos" begin
    ##= Download Files =##
    paths_kernel = [
        "fk/hera_v08.tf",
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
    et_begin = SPICE.utc2et("2027-02-18T00:00:00")
    et_end   = SPICE.utc2et("2027-02-18T01:00:00")
    step     = 300
    et_range = et_begin : step : et_end
    @show et_range
    @show length(et_range)

    # Indices of et_range to be saved.
    # Save only the last rotation.
    save_range = findall(et_range .> et_range[end] - 7.63262 * 3600)
    @show save_range[begin]
    @show save_range[end]
    @show length(save_range)

    # Position 
    sun_d1 = [SVector{3}(SPICE.spkpos("SUN", et, "DIDYMOS_FIXED", "None", "DIDYMOS")[1])*1000 for et in et_range]
    sun_d2 = [SVector{3}(SPICE.spkpos("SUN", et, "DIMORPHOS_FIXED", "None", "DIMORPHOS")[1])*1000 for et in et_range]
    d1_d2 = [SVector{3}(SPICE.spkpos("DIDYMOS", et, "DIMORPHOS_FIXED", "None", "DIMORPHOS")[1])*1000 for et in et_range]
    d2_d1 = [SVector{3}(SPICE.spkpos("DIMORPHOS", et, "DIDYMOS_FIXED", "None", "DIDYMOS")[1])*1000 for et in et_range]
    # Transformation matrix
    D1_TO_D2 = [RotMatrix{3}(SPICE.pxform("DIDYMOS_FIXED", "DIMORPHOS_FIXED", et)) for et in et_range]
    D2_TO_D1 = [RotMatrix{3}(SPICE.pxform("DIMORPHOS_FIXED", "DIDYMOS_FIXED", et)) for et in et_range]
    D1_TO_J2000 = [RotMatrix{3}(SPICE.pxform("DIDYMOS_FIXED", "J2000", et)) for et in et_range]
    D2_TO_J2000 = [RotMatrix{3}(SPICE.pxform("DIMORPHOS_FIXED", "J2000", et)) for et in et_range]
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

    ##= TPM =##
    thermo_params = AsteroidThermoPhysicalModels.thermoparams(  # [Michel+2016; Naidu+2020]
        A_B     = 0.059,  # Bolometric Bond albedo
        A_TH    = 0.0,
        k       = 0.125,
        ρ       = 2170.,
        Cp      = 600.,
        ε       = 0.9,
        t_begin = et_range[begin],
        t_end   = et_range[end],
        Nt      = length(et_range),
        z_max   = 0.6,
        Nz      = 41,
        P       = SPICE.convrt(AsteroidThermoPhysicalModels.DIDYMOS[:P], "hours", "seconds"),
    )

    AsteroidThermoPhysicalModels.init_temperature_zero!(shape1, thermo_params)
    AsteroidThermoPhysicalModels.init_temperature_zero!(shape2, thermo_params)

    # Run TPM and save the result
    savepath = joinpath("TPM_Didymos.jld2")
    shapes = (shape1, shape2)
    suns = (sun_d1, sun_d2)
    AsteroidThermoPhysicalModels.run_TPM!(shapes, et_range, suns, D2_TO_D1, d2_d1, thermo_params, savepath, [:surf_temps, :forces, :torques])
end
