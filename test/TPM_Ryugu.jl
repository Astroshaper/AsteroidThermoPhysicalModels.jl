# See https://github.com/Astroshaper/Astroshaper-examples/tree/main/TPM_Ryugu for more information.
@testset "TPM_Ryugu" begin
    msg = """\n
    ┏━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━┓
    |                    Test: TPM_Ryugu                     |
    ┗━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━┛
    """
    println(msg)

    ##= Download Files =##
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
        url_kernel = "https://data.darts.isas.jaxa.jp/pub/hayabusa2/spice_bundle/spice_kernels/$(path_kernel)"
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

    ##= Load data with SPICE =##
    for path_kernel in paths_kernel
        filepath = joinpath("kernel", path_kernel)
        SPICE.furnsh(filepath)
    end

    ##= Ephemerides =##
    et_begin = SPICE.utc2et("2018-07-01T00:00:00")
    et_end   = SPICE.utc2et("2018-07-01T01:00:00")
    step     = 76.3262  # Rotation of 1 deg
    et_range = et_begin : step : et_end
    @show et_range
    @show length(et_range)

    """
    - `time` : Ephemeris times
    - `sun`  : Sun's position in the RYUGU_FIXED frame
    """
    ephem = (
        time = collect(et_range),
        sun  = [SVector{3}(SPICE.spkpos("SUN", et, "RYUGU_FIXED", "None", "RYUGU")[1]) * 1000 for et in et_range],
    )

    SPICE.kclear()

    ##= Load obj file =##
    path_obj = joinpath("shape", "ryugu_test.obj")   # Small model for test
    path_jld = joinpath("shape", "ryugu_test.jld2")  # Small model for test
    # path_obj = joinpath("shape", "SHAPE_SFM_49k_v20180804.obj")
    # path_jld = joinpath("shape", "SHAPE_SFM_49k_v20180804.jld2")
    if isfile(path_jld) && ENABLE_JLD
        shape = AsteroidThermoPhysicalModels.load_shape_jld(path_jld)
    else
        shape = AsteroidThermoPhysicalModels.load_shape_obj(path_obj; scale=1000, find_visible_facets=true)
        AsteroidThermoPhysicalModels.save_shape_jld(path_jld, shape)
    end

    ##= Thermal properties =##
    P = SPICE.convrt(7.63262, "hours", "seconds")
    k = 0.1
    ρ = 1270.0
    Cₚ = 600.0
    
    l = AsteroidThermoPhysicalModels.thermal_skin_depth(P, k, ρ, Cₚ)
    Γ = AsteroidThermoPhysicalModels.thermal_inertia(k, ρ, Cₚ)

    thermo_params = AsteroidThermoPhysicalModels.thermoparams(
        P       = P,
        l       = l,
        Γ       = Γ,
        A_B     = 0.04,  # Bolometric Bond albedo
        A_TH    = 0.0,
        ε       = 1.0,
        t_begin = et_range[begin],
        t_end   = et_range[end],
        Nt      = length(et_range),
        z_max   = 0.6,
        Nz      = 41,
    )

    println(thermo_params)

    # Run TPM and save the result
    AsteroidThermoPhysicalModels.init_temperature!(shape, thermo_params, 200.)
    savepath = joinpath("TPM_Ryugu.jld2")
    AsteroidThermoPhysicalModels.run_TPM!(shape, thermo_params, ephem, savepath)
end
