
@testset "heat_conduction_1D" begin
    msg = """

    ⋅----------------------------------------⋅
    |        Test: heat_conduction_1D        |
    ⋅----------------------------------------⋅
    """
    println(msg)

    ##= Shape model =##
    path_obj = joinpath("shape", "single_face.obj")
    shape = AsteroidThermoPhysicalModels.load_shape_obj(path_obj)

    ##= Seeting of time step =##
    et_begin = 0.0
    et_end   = 1.0
    step     = 0.001
    et_range = et_begin : step : et_end

    ##= Thermal properties =##
    P  = 1.0
    k  = 1.0
    ρ  = 1.0
    Cₚ = 1.0
    
    l = AsteroidThermoPhysicalModels.thermal_skin_depth(P, k, ρ, Cₚ)
    Γ = AsteroidThermoPhysicalModels.thermal_inertia(k, ρ, Cₚ)

    thermo_params = AsteroidThermoPhysicalModels.thermoparams(
        P       = P,
        l       = l,
        Γ       = Γ,
        A_B     = 0.0,
        A_TH    = 0.0,
        ε       = 1.0,
        t_begin = et_range[begin],
        t_end   = et_range[end],
        Nt      = length(et_range),
        z_max   = 1.0,
        Nz      = 11,
    )

    println(thermo_params)

    ##= Initial temperature =##
    AsteroidThermoPhysicalModels.init_temperature_zero!(shape, thermo_params)

    for n in 1:size(shape.temperature, 1)
        z = thermo_params.Δz * (n - 1)
        if  z < 0.5
            shape.temperature[n, :, :] .= 2z
        else
            shape.temperature[n, :, :] .= 2(1 - z)
        end
    end

    for nₜ in eachindex(et_range)
        et = et_range[nₜ]

        nₜ == length(et_range) && break  # Stop to update the temperature at the final step
        
        λ = thermo_params.λ
        Tⱼ   = @views shape.temperature[:, :, nₜ  ]
        Tⱼ₊₁ = @views shape.temperature[:, :, nₜ+1]
        
        ## Forward Euler method
        @. Tⱼ₊₁[begin+1:end-1, :] = @views (1-2λ')*Tⱼ[begin+1:end-1, :] + λ'*(Tⱼ[begin+2:end, :] + Tⱼ[begin:end-2, :])
        
        ## Isothermal boundary condition
        AsteroidThermoPhysicalModels.update_surface_temperature!(shape, thermo_params, nₜ+1, AsteroidThermoPhysicalModels.Isothermal)
        AsteroidThermoPhysicalModels.update_bottom_temperature!(shape, nₜ+1, AsteroidThermoPhysicalModels.Isothermal)
    end

    savepath = "heat_conduction_1D.jld2"
    jldsave(savepath; shape, et_range, thermo_params)
end
