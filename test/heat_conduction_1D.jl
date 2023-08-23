
@testset "heat_conduction_1D" begin
    msg = """\n
    ┏━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━┓
    |                Test: heat_conduction_1D                |
    ┗━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━┛
    """
    println(msg)

    ##= Shape model =##
    path_obj = joinpath("shape", "single_face.obj")
    shape = AsteroidThermoPhysicalModels.load_shape_obj(path_obj)

    ##= Seeting of time step =##
    et_begin = 0.0
    et_end   = 1.0
    step     = 0.4e-4
    et_range = et_begin : step : et_end

    ephem = (
        time = collect(et_range),
    )

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
        Nz      = 101,
    )

    println(thermo_params)

    ##= TPMs with different solvers =##
    stpm_FE = AsteroidThermoPhysicalModels.SingleTPM(shape, thermo_params; SELF_SHADOWING=false, SELF_HEATING=false, SOLVER=AsteroidThermoPhysicalModels.ForwardEulerSolver())
    stpm_BE = AsteroidThermoPhysicalModels.SingleTPM(shape, thermo_params; SELF_SHADOWING=false, SELF_HEATING=false, SOLVER=AsteroidThermoPhysicalModels.BackwardEulerSolver(thermo_params))
    stpm_CN = AsteroidThermoPhysicalModels.SingleTPM(shape, thermo_params; SELF_SHADOWING=false, SELF_HEATING=false, SOLVER=AsteroidThermoPhysicalModels.CrankNicolsonSolver(thermo_params))

    ##= Initial temperature =##
    T₀(x) = x < 0.5 ? 2x : 2(1 - x)
    xs = [thermo_params.Δz * (nz-1) for nz in 1:thermo_params.Nz]
    Ts = [T₀(x) for x in xs]

    stpm_FE.temperature .= Ts
    stpm_BE.temperature .= Ts
    stpm_CN.temperature .= Ts

    ##= Run TPM =##
    @time for nₜ in eachindex(ephem.time)
        nₜ == length(et_range) && break  # Stop to update the temperature at the final step
        
        AsteroidThermoPhysicalModels.forward_euler!(stpm_FE, nₜ)
        AsteroidThermoPhysicalModels.backward_euler!(stpm_BE, nₜ)
        AsteroidThermoPhysicalModels.crank_nicolson!(stpm_CN, nₜ)
        
        ## Isothermal boundary condition
        # AsteroidThermoPhysicalModels.update_surface_temperature!(stpm, nₜ+1, AsteroidThermoPhysicalModels.Isothermal)
        # AsteroidThermoPhysicalModels.update_bottom_temperature!(stpm, nₜ+1, AsteroidThermoPhysicalModels.Isothermal)

        stpm_FE.temperature[begin, :, nₜ+1] .= 0
        stpm_FE.temperature[end  , :, nₜ+1] .= 0
        stpm_BE.temperature[begin, :, nₜ+1] .= 0
        stpm_BE.temperature[end  , :, nₜ+1] .= 0
        stpm_CN.temperature[begin, :, nₜ+1] .= 0
        stpm_CN.temperature[end  , :, nₜ+1] .= 0
    end

    ##= Save data =##
    df = DataFrames.DataFrame(
        x = xs,
        T_FE_100Δt = stpm_FE.temperature[:, 1, 101],  # t =  4 ms <- 0.4e-4 * 100
        T_BE_100Δt = stpm_BE.temperature[:, 1, 101],  # t =  4 ms
        T_CN_100Δt = stpm_CN.temperature[:, 1, 101],  # t =  4 ms
        T_FE_200Δt = stpm_FE.temperature[:, 1, 201],  # t =  8 ms
        T_BE_200Δt = stpm_BE.temperature[:, 1, 201],  # t =  8 ms
        T_CN_200Δt = stpm_CN.temperature[:, 1, 201],  # t =  8 ms
        T_FE_400Δt = stpm_FE.temperature[:, 1, 401],  # t = 16 ms
        T_BE_400Δt = stpm_BE.temperature[:, 1, 401],  # t = 16 ms
        T_CN_400Δt = stpm_CN.temperature[:, 1, 401],  # t = 16 ms
        T_FE_800Δt = stpm_FE.temperature[:, 1, 801],  # t = 32 ms
        T_BE_800Δt = stpm_BE.temperature[:, 1, 801],  # t = 32 ms
        T_CN_800Δt = stpm_CN.temperature[:, 1, 801],  # t = 32 ms
    )
    jldsave("heat_conduction_1D.jld2"; df)
end
