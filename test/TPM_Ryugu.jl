# See https://github.com/Astroshaper/Astroshaper-examples/tree/main/TPM_Ryugu for more information.
@testset "TPM_Ryugu" begin
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

    ##= Load obj file =##
    path_obj = joinpath("shape", "SHAPE_SFM_49k_v20180804.obj")
    path_jld = joinpath("shape", "SHAPE_SFM_49k_v20180804.jld2")
    if isfile(path_jld)
        shape = AsteroidThermoPhysicalModels.ShapeModel(path_jld)
    else
        shape = AsteroidThermoPhysicalModels.ShapeModel(path_obj; scale=1000, find_visible_facets=true, save_shape=true)
    end

    ##= TPM =##
    thermo_params = AsteroidThermoPhysicalModels.thermoparams(
        A_B     = 0.04,  # Bolometric Bond albedo
        A_TH    = 0.0,
        k       = 0.1,
        ρ       = 1270.0,
        Cp      = 600.0,
        ε       = 1.0,
        t_begin = et_range[begin],
        t_end   = et_range[end],
        Nt      = length(et_range),
        z_max   = 0.6,
        Nz      = 41,
        P       = 7.63262 * 3600,
    )

    # Run TPM and save the result
    savepath = joinpath("TPM_Ryugu.jld2")
    AsteroidThermoPhysicalModels.run_TPM!(shape, et_range, sun_ryugu, thermo_params, savepath, save_range)

    JLD2.jldopen(savepath, "r+") do file
        file["RYUGU_TO_J2000"] = RYUGU_TO_J2000[save_range]
    end
end
