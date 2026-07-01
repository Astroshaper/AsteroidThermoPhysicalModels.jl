#=
heat_conduction_1D.jl

Tests for 1D heat conduction solvers:
- Compares Explicit Euler, Implicit Euler, and Crank-Nicolson methods
- Validates against analytical solutions for isothermal boundaries
- Checks numerical accuracy and stability
- Ensures all three methods converge to the same solution
=#

@testset "heat_conduction_1D" begin
    msg = """\n
    ┏━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━┓
    |                Test: heat_conduction_1D                |
    ┗━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━┛
    """
    println(msg)

    ## --- Shape model ---
    path_obj = joinpath("shape", "single_face.obj")
    shape = load_shape_obj(path_obj)
    n_face = length(shape.faces)  # Number of faces

    ## --- Seeting of time step ---
    et_range = range(0.0, 1.0; step=1e-5)

    times = collect(et_range)

    ## --- Thermal properties ---
    k  = 0.1     # Thermal conductivity [W/m/K]
    ρ  = 1270.0  # Density [kg/m³]
    Cₚ = 600.0   # Heat capacity [J/kg/K]

    R_vis = 0.0  # Reflectance in visible light [-]
    R_ir  = 0.0  # Reflectance in thermal infrared [-]
    ε     = 1.0  # Emissivity [-]

    z_max   = 1.0  # Depth of the lower boundary of a heat conduction equation [m]
    n_depth = 101  # Number of depth steps

    thermo_params = ThermoParams(k, ρ, Cₚ, R_vis, R_ir, ε)
    grid_params   = GridParams(; z_max, n_depth)

    ## --- Build states with different solvers ---
    BC_UPPER = IsothermalBoundaryCondition(0)
    BC_LOWER = IsothermalBoundaryCondition(0)

    problem = SingleAsteroidThermoPhysicalProblem(shape, thermo_params, grid_params;
        with_self_shadowing      = false,
        with_self_heating        = false,
        upper_boundary_condition = BC_UPPER,
        lower_boundary_condition = BC_LOWER,
    )

    state_EE = AsteroidThermoPhysicalModels._build_single_state(problem, ExplicitEuler())
    state_IE = AsteroidThermoPhysicalModels._build_single_state(problem, ImplicitEuler())
    state_CN = AsteroidThermoPhysicalModels._build_single_state(problem, CrankNicolson())

    ## --- Initial temperature ---
    T₀(x) = x < 0.5 ? 2x : 2(1 - x)
    xs = [grid_params.Δz * (i-1) for i in 1:grid_params.n_depth]
    Ts = [T₀(x) for x in xs]

    state_EE.temperature .= Ts
    state_IE.temperature .= Ts
    state_CN.temperature .= Ts

    ## --- Run heat conduction solvers ---
    for i_time in eachindex(times)
        i_time == length(et_range) && break  # Stop to update the temperature at the final step
        Δt = times[i_time+1] - times[i_time]

        AsteroidThermoPhysicalModels.explicit_euler!(state_EE, Δt)
        AsteroidThermoPhysicalModels.implicit_euler!(state_IE, Δt)
        AsteroidThermoPhysicalModels.crank_nicolson!(state_CN, Δt)
    end

    """
        relative_error(a, b) -> δ

    Calculate the relative error between `a` and `b`.

    # Arguments
    - `a` : Scholar value
    - `b` : Scholar value

    # Return
    - `δ` : Relative error between `a` and `b`
    """
    function relative_error(a::Real, b::Real)
        eps = 1e-10  # Small values to avoid zero division
        δ = abs(a - b) / (abs(b) + eps)

        return δ
    end

    @testset "Solver comparison" begin
        ## Maximum relative error between solvers
        δ_max_EE_IE = maximum(relative_error.(state_EE.temperature, state_IE.temperature))
        δ_max_EE_CN = maximum(relative_error.(state_EE.temperature, state_CN.temperature))
        δ_max_IE_CN = maximum(relative_error.(state_IE.temperature, state_CN.temperature))

        ## Allowable relative error
        @test δ_max_EE_IE < 0.01
        @test δ_max_EE_CN < 0.01
        @test δ_max_IE_CN < 0.01

        println("Maximum relative errors:")
        println("    - Explicit Euler vs. Implicit Euler : ", δ_max_EE_IE)
        println("    - Explicit Euler vs. Crank-Nicolson : ", δ_max_EE_CN)
        println("    - Implicit Euler vs. Crank-Nicolson : ", δ_max_IE_CN)
    end

    @testset "Numercal vs. analytical solution (isothermal boundary condition)" begin
        T_analytical = similar(Ts)
        α = AsteroidThermoPhysicalModels.thermal_diffusivity(k, ρ, Cₚ)  # Thermal diffusivity [m²/s]
        for (i, x) in enumerate(xs)
            T_analytical[i] = AsteroidThermoPhysicalModels.analytical_solution_isothermal(x, times[end], z_max, α; n_max=100)
        end

        # Maximum relative error between numerical and analytical solutions
        δ_max_EE = maximum(relative_error.(state_EE.temperature[:, 1], T_analytical))
        δ_max_IE = maximum(relative_error.(state_IE.temperature[:, 1], T_analytical))
        δ_max_CN = maximum(relative_error.(state_CN.temperature[:, 1], T_analytical))

        ## Allowable relative error
        @test δ_max_EE < 0.01
        @test δ_max_IE < 0.01
        @test δ_max_CN < 0.01

        println("Maximum relative errors:")
        println("    - Explicit Euler vs. Isothermal analytical solution: ", δ_max_EE)
        println("    - Implicit Euler vs. Isothermal analytical solution: ", δ_max_IE)
        println("    - Crank-Nicolson vs. Isothermal analytical solution: ", δ_max_CN)
    end
end
