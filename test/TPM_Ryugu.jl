# See https://github.com/Astroshaper/Astroshaper-examples/tree/main/TPM_Ryugu for more information.
@testset "TPM_Ryugu" begin
    ##= Download Files =##
    path_kernels = [
        "spice_kernels/lsk/naif0012.tls",
        "spice_kernels/pck/hyb2_ryugu_shape_v20190328.tpc",
        "spice_kernels/fk/hyb2_ryugu_v01.tf",
        "spice_kernels/spk/2162173_Ryugu.bsp",
    ]

    for path_kernel in path_kernels
        mkpath(dirname(path_kernel))
        url_kernel = joinpath("https://data.darts.isas.jaxa.jp/pub/hayabusa2/spice_bundle", path_kernel)
        isfile(path_kernel) || Downloads.download(url_kernel, path_kernel)
    end

    path_obj = "SHAPE_SFM_49k_v20180804.obj"
    if !isfile(path_obj)
        url_obj = "https://data.darts.isas.jaxa.jp/pub/hayabusa2/paper/Watanabe_2019/SHAPE_SFM_49k_v20180804.obj"
        Downloads.download(url_obj, path_obj)
    end

    ##= Load data with SPICE =##
    for path_kernel in path_kernels
        SPICE.furnsh(path_kernel)
    end

    et_start = SPICE.utc2et("2018-07-01T00:00:00")
    et_end   = SPICE.utc2et("2018-07-01T01:00:00")
    step     = 76.3262  # Rotation of 1 deg
    et_range = et_start : step : et_end
    @show et_range
    @show length(et_range)
    
    # Indices of et_range to be saved.
    # Save only the last rotation.
    save_range = findall(et_range .> et_range[end] - 7.63262 * 3600)
    @show save_range[begin]
    @show save_range[end]
    @show length(save_range);
    
    # Sun's position in the RYUGU_FIXED frame
    sun_ryugu = [SPICE.spkpos("SUN", et, "RYUGU_FIXED", "None", "RYUGU")[1]*1000 for et in et_range]
    # Transformation matrix from RYUGU_FIXED to J2000
    RYUGU_TO_J2000 = [SPICE.pxform("RYUGU_FIXED", "J2000", et) for et in et_range]
    SPICE.kclear()

    ##= Load obj file =##
    path_jld = splitext(path_obj)[1]*".jld2"
    if isfile(path_jld)
        shape = ThermoPhysicalModeling.ShapeModel(path_jld; scale=1000, find_visible_facets=true, save_shape=true)
    else
        shape = ThermoPhysicalModeling.ShapeModel(path_obj; scale=1000, find_visible_facets=true, save_shape=true)
    end

    ##= TPM =##
    thermo_params = ThermoPhysicalModeling.ThermoParams(
        A_B   = 0.04,  # Bolometric Bond albedo
        A_TH  = 0.0,
        k     = 0.1,
        ρ     = 1270.0,
        Cp    = 600.0,
        ϵ     = 1.0,
        t_bgn = et_range[begin],
        t_end = et_range[end],
        Nt    = length(et_range),
        z_max = 0.6,
        Nz    = 41,
        P     = 7.63262 * 3600,
    )
    # Run TPM and save the result
    savepath = "TPM_Ryugu.jld2"
    ThermoPhysicalModeling.run_TPM!(shape, et_range, sun_ryugu, thermo_params, savepath, save_range)
    JLD2.jldopen(savepath, "r+") do file
        file["RYUGU_TO_J2000"] = RYUGU_TO_J2000[save_range]
    end
end