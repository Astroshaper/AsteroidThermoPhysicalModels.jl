#=
test_boundary_conditions.jl

Tests for various boundary conditions in heat conduction:
- Radiation boundary condition (surface energy balance)
- Insulation boundary condition (zero flux)
- Isothermal boundary condition (fixed temperature)
Validates both upper and lower boundary implementations.
=#

@testset "Boundary Conditions" begin
    msg = """
    ┏━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━┓
    |               Test: Boundary Conditions                |
    ┗━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━┛
    """
    println(msg)
    # Simple 1D heat conduction test with non-zero boundary conditions

    # Shape model (single face)
    path_obj = joinpath(@__DIR__, "shape", "single_face.obj")
    shape = load_shape_obj(path_obj)

    # Time settings
    et_range = range(0.0, 0.1; step=1e-4)
    times = collect(et_range)

    # Thermal properties
    k  = 1.0                    # Thermal conductivity [W/m/K]
    ρ  = 1000.0                 # Density [kg/m³]
    Cₚ = 1000.0                 # Heat capacity [J/kg/K]
    R_vis = 0.0                 # Reflectance in visible light [-]
    R_ir  = 0.0                 # Reflectance in thermal infrared [-]
    ε     = 1.0                 # Emissivity [-]
    z_max = 1.0                 # Depth of the lower boundary of a heat conduction equation [m]
    n_depth = 21                # Number of depth steps
    Δz = z_max / (n_depth - 1)  # Depth step width [m]

    thermo_params = ThermoParams(k, ρ, Cₚ, R_vis, R_ir, ε, z_max, Δz, n_depth)

    @testset "Isothermal BC with non-zero temperatures" begin
        # Test with T_upper = 100, T_lower = 50
        BC_UPPER = AsteroidThermoPhysicalModels.IsothermalBoundaryCondition(100.0)
        BC_LOWER = AsteroidThermoPhysicalModels.IsothermalBoundaryCondition(50.0)

        problem = SingleAsteroidThermoPhysicalProblem(shape, thermo_params;
            with_self_shadowing      = false,
            with_self_heating        = false,
            upper_boundary_condition = BC_UPPER,
            lower_boundary_condition = BC_LOWER,
        )
        state = AsteroidThermoPhysicalModels._build_single_state(problem, ImplicitEuler())

        # Initial uniform temperature
        AsteroidThermoPhysicalModels.init_temperature!(state, 75.0)

        # Run simulation
        for i_time in 1:(length(times)-1)
            Δt = times[i_time+1] - times[i_time]
            AsteroidThermoPhysicalModels.update_temperature!(state, Δt)
        end

        # Check boundary conditions are satisfied
        @test isapprox(state.temperature[begin, 1], 100.0, rtol=1e-6)
        @test isapprox(state.temperature[end, 1], 50.0, rtol=1e-6)

        # Check monotonic temperature profile (should decrease from top to bottom)
        temps = state.temperature[:, 1]
        @test all(temps[i] >= temps[i+1] for i in 1:(length(temps)-1))
    end

    @testset "Insulation BC at lower boundary" begin
        BC_UPPER = AsteroidThermoPhysicalModels.IsothermalBoundaryCondition(100.0)
        BC_LOWER = AsteroidThermoPhysicalModels.InsulationBoundaryCondition()

        problem = SingleAsteroidThermoPhysicalProblem(shape, thermo_params;
            with_self_shadowing      = false,
            with_self_heating        = false,
            upper_boundary_condition = BC_UPPER,
            lower_boundary_condition = BC_LOWER,
        )
        state = AsteroidThermoPhysicalModels._build_single_state(problem, ImplicitEuler())

        # Initial uniform temperature
        AsteroidThermoPhysicalModels.init_temperature!(state, 50.0)

        # Run simulation
        for i_time in 1:(length(times)-1)
            Δt = times[i_time+1] - times[i_time]
            AsteroidThermoPhysicalModels.update_temperature!(state, Δt)
        end

        # Check upper boundary condition
        @test isapprox(state.temperature[begin, 1], 100.0, rtol=1e-6)

        # Check insulation boundary condition (zero gradient)
        @test isapprox(state.temperature[end, 1], state.temperature[end-1, 1], rtol=1e-3)
    end
end
