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
    shape = AsteroidThermoPhysicalModels.load_shape_obj(path_obj)
    n_face = length(shape.faces)  # Number of faces

    ## --- Seeting of time step ---
    et_range = range(0.0, 1.0; step=1e-5)

    ephem = (
        time = collect(et_range),
    )

    ## --- Thermal properties ---
    k  = 0.1     # Thermal conductivity [W/m/K]
    ρ  = 1270.0  # Density [kg/m³]
    Cₚ = 600.0   # Heat capacity [J/kg/K]

    R_vis = 0.0  # Reflectance in visible light [-]
    R_ir  = 0.0  # Reflectance in thermal infrared [-]
    ε     = 1.0  # Emissivity [-]

    z_max   = 1.0  # Depth of the lower boundary of a heat conduction equation [m]
    n_depth = 101  # Number of depth steps
    Δz = z_max / (n_depth - 1)  # Depth step width [m]

    thermo_params = AsteroidThermoPhysicalModels.ThermoParams(k, ρ, Cₚ, R_vis, R_ir, ε, z_max, Δz, n_depth)

    ## --- TPMs with different solvers ---
    SELF_SHADOWING = false
    SELF_HEATING   = false
    BC_UPPER       = AsteroidThermoPhysicalModels.IsothermalBoundaryCondition(0)
    BC_LOWER       = AsteroidThermoPhysicalModels.IsothermalBoundaryCondition(0)

    stpm_EE = AsteroidThermoPhysicalModels.SingleAsteroidTPM(shape, thermo_params; SELF_SHADOWING, SELF_HEATING, BC_UPPER, BC_LOWER, SOLVER=AsteroidThermoPhysicalModels.ExplicitEulerSolver(thermo_params))
    stpm_IE = AsteroidThermoPhysicalModels.SingleAsteroidTPM(shape, thermo_params; SELF_SHADOWING, SELF_HEATING, BC_UPPER, BC_LOWER, SOLVER=AsteroidThermoPhysicalModels.ImplicitEulerSolver(thermo_params))
    stpm_CN = AsteroidThermoPhysicalModels.SingleAsteroidTPM(shape, thermo_params; SELF_SHADOWING, SELF_HEATING, BC_UPPER, BC_LOWER, SOLVER=AsteroidThermoPhysicalModels.CrankNicolsonSolver(thermo_params))

    ## --- Initial temperature ---
    T₀(x) = x < 0.5 ? 2x : 2(1 - x)
    xs = [thermo_params.Δz * (i-1) for i in 1:thermo_params.n_depth]
    Ts = [T₀(x) for x in xs]

    stpm_EE.temperature .= Ts
    stpm_IE.temperature .= Ts
    stpm_CN.temperature .= Ts

    ## --- Run TPM ---
    for i_time in eachindex(ephem.time)
        i_time == length(et_range) && break  # Stop to update the temperature at the final step
        Δt = ephem.time[i_time+1] - ephem.time[i_time]
        
        AsteroidThermoPhysicalModels.explicit_euler!(stpm_EE, Δt)
        AsteroidThermoPhysicalModels.implicit_euler!(stpm_IE, Δt)
        AsteroidThermoPhysicalModels.crank_nicolson!(stpm_CN, Δt)
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
        δ_max_EE_IE = maximum(relative_error.(stpm_EE.temperature, stpm_IE.temperature))
        δ_max_EE_CN = maximum(relative_error.(stpm_EE.temperature, stpm_CN.temperature))
        δ_max_IE_CN = maximum(relative_error.(stpm_IE.temperature, stpm_CN.temperature))
        
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
            T_analytical[i] = AsteroidThermoPhysicalModels.analytical_solution_isothermal(x, ephem.time[end], z_max, α; n_max=100)
        end

        # Maximum relative error between numerical and analytical solutions
        δ_max_EE = maximum(relative_error.(stpm_EE.temperature[:, 1], T_analytical))
        δ_max_IE = maximum(relative_error.(stpm_IE.temperature[:, 1], T_analytical))
        δ_max_CN = maximum(relative_error.(stpm_CN.temperature[:, 1], T_analytical))

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
