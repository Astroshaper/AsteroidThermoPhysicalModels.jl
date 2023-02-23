# See https://github.com/MasanoriKanamaru/Astroshaper-examples/tree/main/TPM_Didymos for more information.
@testset "TPM_Didymos" begin
    ##= Download Files =##
    path_meta_new = "hera_study_PO_EMA_2024.tm"
    if !isfile(path_meta_new)
        # Dowonload files for SPICE
        run(Git.git(["clone", "https://s2e2.cosmos.esa.int/bitbucket/scm/spice_kernels/hera.git"]))
    
        # Update path
        path_meta_old = "hera/kernels/mk/hera_study_PO_EMA_2024.tm"
        cp(path_meta_old, path_meta_new, force=true)
        script = read(path_meta_new, String)
        script = replace(script, "PATH_VALUES       = ( '..' )"=>"PATH_VALUES     = ('$(abspath(joinpath("hera", "kernels")))')")
        write(path_meta_new, script)
    end

    ##= Load data with SPICE =##
    SPICE.furnsh(path_meta_new)
    et_start = SPICE.utc2et("2027-02-18T00:00:00")
    et_end   = SPICE.utc2et("2027-02-18T01:00:00")
    step     = 300
    et_range = et_start : step : et_end
    @show et_range
    @show length(et_range)
    
    # Indices of et_range to be saved.
    # Save only the last rotation.
    save_range = findall(et_range .> et_range[end] - 7.63262 * 3600)
    @show save_range[begin]
    @show save_range[end]
    @show length(save_range)
    
    # Position 
    sun_d1 = [SPICE.spkpos("SUN", et, "DIDYMOS_FIXED", "None", "DIDYMOS")[1]*1000 for et in et_range]
    sun_d2 = [SPICE.spkpos("SUN", et, "DIMORPHOS_FIXED", "None", "DIMORPHOS")[1]*1000 for et in et_range]
    d1_d2 = [SPICE.spkpos("DIDYMOS", et, "DIMORPHOS_FIXED", "None", "DIMORPHOS")[1]*1000 for et in et_range]
    d2_d1 = [SPICE.spkpos("DIMORPHOS", et, "DIDYMOS_FIXED", "None", "DIDYMOS")[1]*1000 for et in et_range]
    # Transformation matrix
    D1_TO_D2 = [SPICE.pxform("DIDYMOS_FIXED", "DIMORPHOS_FIXED", et) for et in et_range]
    D2_TO_D1 = [SPICE.pxform("DIMORPHOS_FIXED", "DIDYMOS_FIXED", et) for et in et_range]
    D1_TO_J2000 = [SPICE.pxform("DIDYMOS_FIXED", "J2000", et) for et in et_range]
    D2_TO_J2000 = [SPICE.pxform("DIMORPHOS_FIXED", "J2000", et) for et in et_range]
    SPICE.kclear()

    ##= Load obj file =##
    path_shape1_obj = "g_50677mm_rad_obj_dida_0000n00000_v001.obj"
    path_shape2_obj = "g_06650mm_rad_obj_didb_0000n00000_v001.obj"
    path_shape1_jld = "g_50677mm_rad_obj_dida_0000n00000_v001.jld2"
    path_shape2_jld = "g_06650mm_rad_obj_didb_0000n00000_v001.jld2"
    
    cp(joinpath("hera", "kernels", "dsk", "g_50677mm_rad_obj_dida_0000n00000_v001.obj"), path_shape1_obj, force=true)
    cp(joinpath("hera", "kernels", "dsk", "g_06650mm_rad_obj_didb_0000n00000_v001.obj"), path_shape2_obj, force=true)
    
    if isfile(path_shape1_jld)
        shape1 = ThermoPhysicalModeling.ShapeModel(path_shape1_jld; scale=1000, find_visible_facets=true, save_shape=true)
    else
        shape1 = ThermoPhysicalModeling.ShapeModel(path_shape1_obj; scale=1000, find_visible_facets=true, save_shape=true)
    end
    if isfile(path_shape2_jld)
        shape2 = ThermoPhysicalModeling.ShapeModel(path_shape2_jld; scale=1000, find_visible_facets=true, save_shape=true)
    else
        shape2 = ThermoPhysicalModeling.ShapeModel(path_shape2_obj; scale=1000, find_visible_facets=true, save_shape=true)
    end

    ##= TPM =##
    thermo_params = ThermoPhysicalModeling.ThermoParams(  # [Michel+2016; Naidu+2020]
        A_B   = 0.059,  # Bolometric Bond albedo
        A_TH  = 0.0,
        k     = 0.125,
        ρ     = 2170.,
        Cp    = 600.,
        ϵ     = 0.9,
        t_bgn = et_range[begin],
        t_end = et_range[end],
        Nt    = length(et_range),
        z_max = 0.6,
        Nz    = 41,
        P     = SPICE.convrt(ThermoPhysicalModeling.DIDYMOS[:P], "hours", "seconds"),
    );

    ThermoPhysicalModeling.init_temps_zero!(shape1, thermo_params)
    ThermoPhysicalModeling.init_temps_zero!(shape2, thermo_params)
    
    # Run TPM and save the result
    savepath = "TPM_Didymos.jld2"
    shapes = (shape1, shape2)
    suns = (sun_d1, sun_d2)
    ThermoPhysicalModeling.run_TPM!(shapes, et_range, suns, D2_TO_D1, d2_d1, thermo_params, savepath, [:surf_temps, :forces, :torques])
end
