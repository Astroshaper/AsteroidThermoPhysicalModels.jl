
"""
    analytical_solution_isothermal(x, t, L, α; n_max=100) -> T

Calculate the analytical solution of the 1D heat equation with isothermal boundary conditions.
- Equation: ∂T/∂t = α ∂²T/∂x²
- Domain: 0 ≤ x ≤ L
- Boundary conditions: T(0,t) = T(L,t) = 0
- Initial condition: T₀(x) = x < 0.5L ? 2x/L : 2(1 - x/L)  # Triangular profile
        T₀
        ^ 
      1 |   ・
        |  / \
        | /   \
        |/     \
      0 +---+---+--> x
        0  L/2  L

The solution is given by the Fourier series:
    T(x, t) = Σ Bₙ * sin(nπx/L) * exp(-αn²π²t/L²)

where Bₙ = (2/L) * ∫₀^L T₀(ξ) * sin(nπξ/L) dξ

For the triangular initial condition, the coefficients can be calculated analytically:
    Bₙ = (8/n²π²) * sin(nπ/2)
only for odd `n`. The sum of even-`n` terms is zero due to symmetry.

# Arguments
- `x`     : Position [m]
- `t`     : Time [s]
- `L`     : Length of the domain [m]
- `α`     : Thermal diffusivity [m²/s]
- `n_max` : Number of terms in the Fourier series

# Returns
- `T`     : Temperature [K]
"""
function analytical_solution_isothermal(x, t, L, α; n_max=100)
    T = 0.0
    for n in 1:2:n_max  # Only odd terms contribute for this initial condition
        Bₙ = 8/ (n * π)^2 * sin(n*π/2)
        T += Bₙ * sin(n*π*x/L) * exp(-α*(n*π/L)^2*t)
    end
    return T
end


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
    n_face = length(shape.faces)  # Number of faces

    ##= Seeting of time step =##
    et_range = range(0.0, 1.0; step=0.4e-4)

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

    R_vis = 0.0  # Reflectance in visible light [-]
    R_ir  = 0.0  # Reflectance in thermal infrared [-]
    ε     = 1.0  # Emissivity [-]

    z_max   = 1.0  # Depth of the lower boundary of a heat conduction equation [m]
    n_depth = 101  # Number of depth steps
    Δz = z_max / (n_depth - 1)  # Depth step width [m]

    thermo_params = AsteroidThermoPhysicalModels.ThermoParams(
        P,
        fill(l,     n_face),
        fill(Γ,     n_face),
        fill(R_vis, n_face),
        fill(R_ir,  n_face),
        fill(ε,     n_face),
        z_max,
        Δz,
        n_depth
    )

    ##= TPMs with different solvers =##
    SELF_SHADOWING = false
    SELF_HEATING   = false
    BC_UPPER       = AsteroidThermoPhysicalModels.IsothermalBoundaryCondition(0)
    BC_LOWER       = AsteroidThermoPhysicalModels.IsothermalBoundaryCondition(0)

    stpm_FE = AsteroidThermoPhysicalModels.SingleTPM(shape, thermo_params; SELF_SHADOWING, SELF_HEATING, BC_UPPER, BC_LOWER, SOLVER=AsteroidThermoPhysicalModels.ForwardEulerSolver(thermo_params))
    stpm_BE = AsteroidThermoPhysicalModels.SingleTPM(shape, thermo_params; SELF_SHADOWING, SELF_HEATING, BC_UPPER, BC_LOWER, SOLVER=AsteroidThermoPhysicalModels.BackwardEulerSolver(thermo_params))
    stpm_CN = AsteroidThermoPhysicalModels.SingleTPM(shape, thermo_params; SELF_SHADOWING, SELF_HEATING, BC_UPPER, BC_LOWER, SOLVER=AsteroidThermoPhysicalModels.CrankNicolsonSolver(thermo_params))

    ##= Initial temperature =##
    T₀(x) = x < 0.5 ? 2x : 2(1 - x)
    xs = [thermo_params.Δz * (i-1) for i in 1:thermo_params.n_depth]
    Ts = [T₀(x) for x in xs]

    stpm_FE.temperature .= Ts
    stpm_BE.temperature .= Ts
    stpm_CN.temperature .= Ts

    ##= Run TPM =##
    for i_time in eachindex(ephem.time)
        i_time == length(et_range) && break  # Stop to update the temperature at the final step
        Δt = ephem.time[i_time+1] - ephem.time[i_time]
        
        AsteroidThermoPhysicalModels.forward_euler!(stpm_FE, Δt)
        AsteroidThermoPhysicalModels.backward_euler!(stpm_BE, Δt)
        AsteroidThermoPhysicalModels.crank_nicolson!(stpm_CN, Δt)
    end

    ##= Save data =##
    # df = DataFrames.DataFrame(
    #     x = xs,
        # T_FE_100Δt = stpm_FE.temperature[:, 1, 101],  # t =  4 ms <- 0.4e-4 * 100
        # T_BE_100Δt = stpm_BE.temperature[:, 1, 101],  # t =  4 ms
        # T_CN_100Δt = stpm_CN.temperature[:, 1, 101],  # t =  4 ms
        # T_FE_200Δt = stpm_FE.temperature[:, 1, 201],  # t =  8 ms
        # T_BE_200Δt = stpm_BE.temperature[:, 1, 201],  # t =  8 ms
        # T_CN_200Δt = stpm_CN.temperature[:, 1, 201],  # t =  8 ms
        # T_FE_400Δt = stpm_FE.temperature[:, 1, 401],  # t = 16 ms
        # T_BE_400Δt = stpm_BE.temperature[:, 1, 401],  # t = 16 ms
        # T_CN_400Δt = stpm_CN.temperature[:, 1, 401],  # t = 16 ms
        # T_FE_800Δt = stpm_FE.temperature[:, 1, 801],  # t = 32 ms
        # T_BE_800Δt = stpm_BE.temperature[:, 1, 801],  # t = 32 ms
        # T_CN_800Δt = stpm_CN.temperature[:, 1, 801],  # t = 32 ms
    # )
    # jldsave("heat_conduction_1D.jld2"; df)
end
