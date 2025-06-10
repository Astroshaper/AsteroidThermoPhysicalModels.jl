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
    shape = AsteroidThermoPhysicalModels.load_shape_obj(path_obj)
    
    # Time settings
    et_range = range(0.0, 0.1; step=1e-4)
    ephem = (time = collect(et_range),)
    
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
    
    thermo_params = AsteroidThermoPhysicalModels.ThermoParams(k, ρ, Cₚ, R_vis, R_ir, ε, z_max, Δz, n_depth)
    
    @testset "Isothermal BC with non-zero temperatures" begin
        # Test with T_upper = 100, T_lower = 50
        BC_UPPER = AsteroidThermoPhysicalModels.IsothermalBoundaryCondition(100.0)
        BC_LOWER = AsteroidThermoPhysicalModels.IsothermalBoundaryCondition(50.0)
        
        stpm = AsteroidThermoPhysicalModels.SingleAsteroidTPM(
            shape, thermo_params;
            SELF_SHADOWING=false, 
            SELF_HEATING=false,
            BC_UPPER=BC_UPPER,
            BC_LOWER=BC_LOWER,
            SOLVER=AsteroidThermoPhysicalModels.ImplicitEulerSolver(thermo_params)
        )
        
        # Initial uniform temperature
        AsteroidThermoPhysicalModels.init_temperature!(stpm, 75.0)
        
        # Run simulation
        for i_time in 1:(length(ephem.time)-1)
            Δt = ephem.time[i_time+1] - ephem.time[i_time]
            AsteroidThermoPhysicalModels.update_temperature!(stpm, Δt)
        end
        
        # Check boundary conditions are satisfied
        @test isapprox(stpm.temperature[begin, 1], 100.0, rtol=1e-6)
        @test isapprox(stpm.temperature[end, 1], 50.0, rtol=1e-6)
        
        # Check monotonic temperature profile (should decrease from top to bottom)
        temps = stpm.temperature[:, 1]
        @test all(temps[i] >= temps[i+1] for i in 1:(length(temps)-1))
    end
    
    @testset "Insulation BC at lower boundary" begin
        BC_UPPER = AsteroidThermoPhysicalModels.IsothermalBoundaryCondition(100.0)
        BC_LOWER = AsteroidThermoPhysicalModels.InsulationBoundaryCondition()
        
        stpm = AsteroidThermoPhysicalModels.SingleAsteroidTPM(
            shape, thermo_params;
            SELF_SHADOWING=false, 
            SELF_HEATING=false,
            BC_UPPER=BC_UPPER,
            BC_LOWER=BC_LOWER,
            SOLVER=AsteroidThermoPhysicalModels.ImplicitEulerSolver(thermo_params)
        )
        
        # Initial uniform temperature
        AsteroidThermoPhysicalModels.init_temperature!(stpm, 50.0)
        
        # Run simulation
        for i_time in 1:(length(ephem.time)-1)
            Δt = ephem.time[i_time+1] - ephem.time[i_time]
            AsteroidThermoPhysicalModels.update_temperature!(stpm, Δt)
        end
        
        # Check upper boundary condition
        @test isapprox(stpm.temperature[begin, 1], 100.0, rtol=1e-6)
        
        # Check insulation boundary condition (zero gradient)
        @test isapprox(stpm.temperature[end, 1], stpm.temperature[end-1, 1], rtol=1e-3)
    end
end
