# The following tests are almost the same as `TPM_Ryugu.jl`.
# The only difference is that the thermophysical properties vary depending on the location of the asteroid.
@testset "non-uniform_thermoparams" begin
    msg = """\n
    ┏━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━┓
    |             Test: non-uniform_thermoparams             |
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
    et_begin = SPICE.utc2et("2018-07-01T00:00:00")
    et_end   = SPICE.utc2et("2018-07-01T01:00:00")
    step     = 76.3262  # Rotation of 1 deg
    et_range = et_begin : step : et_end
    @show et_range
    @show length(et_range)

    # Indices of et_range to be saved.
    # Save only the last rotation.
    save_range = findall(et_range .> et_range[end] - 7.63262 * 3600)
    @show save_range[begin]
    @show save_range[end]
    @show length(save_range)

    # Sun's position in the RYUGU_FIXED frame
    sun_ryugu = [SPICE.spkpos("SUN", et, "RYUGU_FIXED", "None", "RYUGU")[1]*1000 for et in et_range]
    # Transformation matrix from RYUGU_FIXED to J2000
    RYUGU_TO_J2000 = [SPICE.pxform("RYUGU_FIXED", "J2000", et) for et in et_range]
    SPICE.kclear()

    ##= Ephemerides =##
    ephem = (
        time = et_range,
        sun  = sun_ryugu,
    )

    ##= Load obj file =##
    path_obj = joinpath("shape", "SHAPE_SFM_49k_v20180804.obj")
    path_jld = joinpath("shape", "SHAPE_SFM_49k_v20180804.jld2")
    if isfile(path_jld) && ENABLE_JLD
        shape = AsteroidThermoPhysicalModels.load_shape_jld(path_jld)
    else
        shape = AsteroidThermoPhysicalModels.load_shape_obj(path_obj; scale=1000, find_visible_facets=true)
        AsteroidThermoPhysicalModels.save_shape_jld(path_jld, shape)
    end

    ##= Thermal properties =##
    """
    When thermophysical properties vary from face to face

    - "Northern" hemisphere:
        - Bond albedo          : A_B = 0.04 [-]
        - Thermal conductivity : k   = 0.1  [W/m/K]
        - Emissivity           : ε   = 1.0  [-]
    - "Southern" hemisphere:
        - Bond albedo          : A_B = 0.1  [-]
        - Thermal conductivity : k   = 0.3  [W/m/K]
        - Emissivity           : ε   = 0.9  [-]
    """

    P  = SPICE.convrt(7.63262, "hours", "seconds")
    k  = [r[3] > 0 ? 0.1 : 0.3  for r in shape.face_centers]
    ρ  = 1270.0
    Cₚ = 600.0
    
    l = AsteroidThermoPhysicalModels.thermal_skin_depth(P, k, ρ, Cₚ)
    Γ = AsteroidThermoPhysicalModels.thermal_inertia(k, ρ, Cₚ)

    thermo_params = AsteroidThermoPhysicalModels.thermoparams(
        P       = P,
        l       = l,
        Γ       = Γ,
        A_B     = [r[3] > 0 ? 0.04 : 0.1 for r in shape.face_centers],
        A_TH    = 0.0,
        ε       = [r[3] > 0 ? 1.0 : 0.9  for r in shape.face_centers],
        t_begin = et_range[begin],
        t_end   = et_range[end],
        Nt      = length(et_range),
        z_max   = 0.6,
        Nz      = 41,
    )

    ##= Run TPM and save the result =##
    AsteroidThermoPhysicalModels.init_temperature!(shape, thermo_params, 200.)
    savepath = joinpath("non-uniform_thermoparams.jld2")
    AsteroidThermoPhysicalModels.run_TPM!(shape, thermo_params, ephem, savepath, save_range)
end
