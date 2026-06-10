#=
heat_conduction.jl

1D heat conduction solvers for thermophysical modeling.
This file implements various numerical methods for solving the heat diffusion equation
in the subsurface of asteroids, including:
- Explicit Euler method (conditionally stable)
- Implicit Euler method (unconditionally stable)
- Crank-Nicolson method (second-order accurate)
- Special handling for zero thermal conductivity
- Various boundary conditions (radiation, insulation, isothermal)
=#

# ╔═══════════════════════════════════════════════════════════════════╗
# ║                      1D heat conduction                           ║
# ╚═══════════════════════════════════════════════════════════════════╝

"""
    update_temperature_zero_conductivity!(state::SingleAsteroidThermoPhysicalState)

Update surface temperature for the zero thermal conductivity case.
When thermal conductivity is zero, there is no heat conduction into the subsurface,
and the surface temperature is determined solely by instantaneous radiative equilibrium.

# Arguments
- `state::SingleAsteroidThermoPhysicalState` : Thermophysical simulation state for a single asteroid

# Mathematical Formula
For each face, the surface temperature T is calculated from:
```
εσT⁴ = (1-Rᵥᵢₛ)F_sun + (1-Rᵥᵢₛ)F_scat + (1-Rᵢᵣ)F_rad
```
where:
- ε : Emissivity
- σ : Stefan-Boltzmann constant
- Rᵥᵢₛ : Reflectance in visible light
- Rᵢᵣ : Reflectance in thermal infrared
- F_sun : Direct solar flux
- F_scat : Scattered light flux
- F_rad : Thermal radiation flux from surrounding surfaces

# Notes
- This function is called when `state.problem.thermo_params.thermal_conductivity` is zero
- The temperature instantly adjusts to balance incoming and outgoing radiation
- No subsurface temperatures are updated (only surface layer)
"""
function update_temperature_zero_conductivity!(state::SingleAsteroidThermoPhysicalState)
    for i_face in eachindex(state.problem.shape.faces)
        R_vis = state.problem.thermo_params.reflectance_vis[i_face]
        R_ir  = state.problem.thermo_params.reflectance_ir[i_face]
        ε     = state.problem.thermo_params.emissivity[i_face]
        εσ    = ε * σ_SB

        F_sun = state.flux_sun[i_face]
        F_scat = state.flux_scat[i_face]
        F_rad = state.flux_rad[i_face]
        F_abs = absorbed_energy_flux(R_vis, R_ir, F_sun, F_scat, F_rad)

        state.temperature[begin, i_face] = (F_abs / εσ)^(1/4)
    end
end

"""
    update_temperature!(state::SingleAsteroidThermoPhysicalState, Δt)

Update the temperature distribution for the next time step by solving the 1D heat conduction equation.
The solver method is determined by `state.solver_cache`, and special handling is applied for zero conductivity.

# Arguments
- `state::SingleAsteroidThermoPhysicalState` : Thermophysical simulation state for a single asteroid
- `Δt::Real` : Time step [s]

# Solver Selection
The function automatically selects the appropriate solver based on `state.solver_cache`:
- `ExplicitEulerCache`: Forward Euler method (conditionally stable, requires λ < 0.5)
- `ImplicitEulerCache`: Backward Euler method (unconditionally stable)
- `CrankNicolsonCache`: Crank-Nicolson method (unconditionally stable, second-order accurate)

# Special Cases
- If thermal conductivity is zero, calls `update_temperature_zero_conductivity!` instead
- The zero-conductivity case uses instantaneous radiative equilibrium

# Mathematical Background
Solves the 1D heat conduction equation:
```
∂T/∂t = α ∂²T/∂z²
```
where α = k/(ρCₚ) is the thermal diffusivity.

# See Also
- `explicit_euler!`, `implicit_euler!`, `crank_nicolson!` for specific solver implementations
- `update_temperature_zero_conductivity!` for the zero-conductivity case
"""
function update_temperature!(state::SingleAsteroidThermoPhysicalState, Δt)
    # Handle zero-conductivity case
    if iszero(state.problem.thermo_params.thermal_conductivity)
        update_temperature_zero_conductivity!(state)
        return
    end
    
    # Non-zero conductivity: use selected solver
    if state.solver_cache isa ExplicitEulerCache
        explicit_euler!(state, Δt)
    elseif state.solver_cache isa ImplicitEulerCache
        implicit_euler!(state, Δt)
    elseif state.solver_cache isa CrankNicolsonCache
        crank_nicolson!(state, Δt)
    else
        error("Unknown solver type: $(typeof(state.solver_cache)). Expected ExplicitEulerCache, ImplicitEulerCache, or CrankNicolsonCache.")
    end
end


"""
    update_temperature!(state::BinaryAsteroidThermoPhysicalState, Δt)

Calculate the temperature for the next time step based on 1D heat conductivity equation.

# Arguments
- `state` : Thermophysical simulation state for a binary asteroid
- `Δt`    : Time step [s]
"""
function update_temperature!(state::BinaryAsteroidThermoPhysicalState, Δt)
    update_temperature!(state.primary, Δt)
    update_temperature!(state.secondary, Δt)
end

# ╔═══════════════════════════════════════════════════════════════════╗
# ║               Solvers of a heat conduction equation               ║
# ╚═══════════════════════════════════════════════════════════════════╝


"""
    explicit_euler!(state::SingleAsteroidThermoPhysicalState, Δt)

Solve the 1D heat conduction equation using the explicit (forward) Euler method.
This method is conditionally stable and requires careful time step selection.

# Arguments
- `state::SingleAsteroidThermoPhysicalState` : Thermophysical simulation state for a single asteroid
- `Δt::Real` : Time step [s]

# Method Properties
- **Time discretization**: Explicit (forward difference)
- **Accuracy**: First-order in time, second-order in space
- **Stability**: Conditionally stable, requires λ = αΔt/Δz² < 0.5

# Discretization
The heat conduction equation ∂T/∂t = α∂²T/∂z² is discretized as:
```
T[i,n+1] = T[i,n] + λ(T[i+1,n] - 2T[i,n] + T[i-1,n])
```
where:
- λ = αΔt/Δz² is the dimensionless time step
- α = k/(ρCₚ) is the thermal diffusivity
- n is the time index, i is the depth index

# Stability Criterion
The method is stable only when λ < 0.5. If this condition is violated, an error is thrown.

# Boundary Conditions
- Upper boundary: Determined by `update_upper_temperature!`
- Lower boundary: Determined by `update_lower_temperature!`

# Performance Notes
- This method is simple and fast but requires small time steps for stability
- Consider using implicit methods for larger time steps

# Errors
- Throws an error if λ ≥ 0.5 (stability violation)
"""
function explicit_euler!(state::SingleAsteroidThermoPhysicalState, Δt)
    T = state.temperature
    n_depth = size(T, 1)
    n_face = size(T, 2)

    for i_face in 1:n_face
        k  = state.problem.thermo_params.thermal_conductivity[i_face]
        ρ  = state.problem.thermo_params.density[i_face]
        Cₚ = state.problem.thermo_params.heat_capacity[i_face]
        Δz = state.problem.thermo_params.Δz

        α = thermal_diffusivity(k, ρ, Cₚ)  # Thermal diffusivity [m²/s]
        λ = α * Δt / Δz^2
        λ ≥ 0.5 && error("Explicit Euler method is unstable because λ = αΔt/Δz² = $λ (must be < 0.5). Consider reducing time step from Δt = $Δt or using an implicit solver.")

        for i_depth in 2:(n_depth-1)
            state.solver_cache.x[i_depth] = (1-2λ)*T[i_depth, i_face] + λ*(T[i_depth+1, i_face] + T[i_depth-1, i_face])  # Predict temperature at next time step
        end

        ## Apply boundary conditions
        update_upper_temperature!(state, i_face)
        update_lower_temperature!(state)

        T[:, i_face] .= state.solver_cache.x  # Copy temperature at next time step
    end
end


"""
    implicit_euler!(state::SingleAsteroidThermoPhysicalState, Δt)

Solve the 1D heat conduction equation using the implicit (backward) Euler method.
This method is unconditionally stable, allowing for larger time steps than explicit methods.

# Arguments
- `state::SingleAsteroidThermoPhysicalState` : Thermophysical simulation state for a single asteroid
- `Δt::Real` : Time step [s]

# Method Properties
- **Time discretization**: Implicit (backward difference)
- **Accuracy**: First-order in time, second-order in space
- **Stability**: Unconditionally stable for any time step size

# Discretization
The heat conduction equation ∂T/∂t = α∂²T/∂z² is discretized as:
```
T[i,n+1] - T[i,n] = λ(T[i+1,n+1] - 2T[i,n+1] + T[i-1,n+1])
```
This leads to a tridiagonal system:
```
-λT[i-1,n+1] + (1+2λ)T[i,n+1] - λT[i+1,n+1] = T[i,n]
```
where λ = αΔt/Δz²

# Solution Method
The resulting tridiagonal system is solved using the Thomas algorithm
(tridiagonal matrix algorithm) for each face.

# Boundary Conditions
Different boundary conditions modify the tridiagonal matrix:
- **Radiation BC**: Special treatment after solving the system
- **Insulation BC**: Modified coefficients at boundaries
- **Isothermal BC**: Direct temperature assignment

# Advantages
- Unconditionally stable - no restriction on time step size
- Allows for larger time steps compared to explicit methods
- More computationally intensive per step but often faster overall

# See Also
- `tridiagonal_matrix_algorithm!` for the solution algorithm
- `update_upper_temperature!`, `update_lower_temperature!` for boundary conditions
"""
function implicit_euler!(state::SingleAsteroidThermoPhysicalState, Δt)
    T = state.temperature
    n_depth = size(T, 1)
    n_face = size(T, 2)

    for i_face in 1:n_face
        k  = state.problem.thermo_params.thermal_conductivity[i_face]
        ρ  = state.problem.thermo_params.density[i_face]
        Cₚ = state.problem.thermo_params.heat_capacity[i_face]
        Δz = state.problem.thermo_params.Δz

        α = thermal_diffusivity(k, ρ, Cₚ)  # Thermal diffusivity [m²/s]
        λ = α * Δt / Δz^2

        # Initialize the tridiagonal matrix coefficients
        state.solver_cache.a .= -λ
        state.solver_cache.a[begin] = 0
        state.solver_cache.a[end]   = 0

        state.solver_cache.b .= 1 + 2λ
        state.solver_cache.b[begin] = 1
        state.solver_cache.b[end]   = 1

        state.solver_cache.c .= -λ
        state.solver_cache.c[begin] = 0
        state.solver_cache.c[end]   = 0

        # Set the right-hand side vector
        state.solver_cache.d .= T[:, i_face]
            
        # Apply lower boundary condition to the matrix system
        if state.problem.lower_boundary_condition isa IsothermalBoundaryCondition
            # T[end] = T_iso
            state.solver_cache.a[end] = 0
            state.solver_cache.b[end] = 1
            state.solver_cache.c[end] = 0
            state.solver_cache.d[end] = state.problem.lower_boundary_condition.T_iso
        elseif state.problem.lower_boundary_condition isa InsulationBoundaryCondition
            # ∂T/∂z = 0 → T[n+1] = T[n]
            state.solver_cache.a[end] = -λ
            state.solver_cache.b[end] = 1 + λ
            state.solver_cache.c[end] = 0
            # d[end] already contains T[end, i_face]
        end

        # Apply upper boundary condition to the matrix system
        if state.problem.upper_boundary_condition isa IsothermalBoundaryCondition
            # T[1] = T_iso
            state.solver_cache.a[begin] = 0
            state.solver_cache.b[begin] = 1
            state.solver_cache.c[begin] = 0
            state.solver_cache.d[begin] = state.problem.upper_boundary_condition.T_iso
        elseif state.problem.upper_boundary_condition isa InsulationBoundaryCondition
            # ∂T/∂z = 0 → T[0] = T[1]
            state.solver_cache.a[begin] = 0
            state.solver_cache.b[begin] = 1 + λ
            state.solver_cache.c[begin] = -λ
            # d[begin] already contains T[begin, i_face]
        elseif state.problem.upper_boundary_condition isa RadiationBoundaryCondition
            # For radiation BC, we keep the standard interior equation
            # and handle it separately after solving
        end
            
        # Solve the tridiagonal system
        tridiagonal_matrix_algorithm!(state)
            
        # Special handling for radiation boundary condition
        if state.problem.upper_boundary_condition isa RadiationBoundaryCondition
            update_upper_temperature!(state, i_face)
        end
            
        # Copy temperature at next time step
        T[:, i_face] .= state.solver_cache.x
    end
end


"""
    crank_nicolson!(state::SingleAsteroidThermoPhysicalState, Δt)

Solve the 1D heat conduction equation using the Crank-Nicolson method.
This method combines the explicit and implicit Euler methods for improved accuracy.

# Arguments
- `state::SingleAsteroidThermoPhysicalState` : Thermophysical simulation state for a single asteroid
- `Δt::Real` : Time step [s]

# Method Properties
- **Time discretization**: Semi-implicit (average of forward and backward differences)
- **Accuracy**: Second-order in both time and space
- **Stability**: Unconditionally stable for any time step size

# Discretization
The heat conduction equation ∂T/∂t = α∂²T/∂z² is discretized using the average
of explicit and implicit schemes:
```
T[i,n+1] - T[i,n] = (α∆t)/(2∆z²) × 
    [(T[i+1,n+1] - 2T[i,n+1] + T[i-1,n+1]) + (T[i+1,n] - 2T[i,n] + T[i-1,n])]
```
This leads to a tridiagonal system:
```
-rT[i-1,n+1] + (1+2r)T[i,n+1] - rT[i+1,n+1] = 
    rT[i-1,n] + (1-2r)T[i,n] + rT[i+1,n]
```
where r = α∆t/(2∆z²)

# Advantages
- Higher accuracy than both explicit and implicit Euler methods
- Unconditionally stable
- Optimal balance between accuracy and computational cost
- Second-order accuracy in both time and space

# Implementation Details
The method requires solving a tridiagonal system at each time step,
similar to the implicit Euler method but with a modified right-hand side
that includes information from the current time step.

# See Also
- `tridiagonal_matrix_algorithm!` for the solution algorithm
- `implicit_euler!`, `explicit_euler!` for comparison with other methods
"""
function crank_nicolson!(state::SingleAsteroidThermoPhysicalState, Δt)
    T = state.temperature
    n_depth = size(T, 1)
    n_face = size(T, 2)

    for i_face in 1:n_face
        k  = state.problem.thermo_params.thermal_conductivity[i_face]
        ρ  = state.problem.thermo_params.density[i_face]
        Cₚ = state.problem.thermo_params.heat_capacity[i_face]
        Δz = state.problem.thermo_params.Δz

        α = thermal_diffusivity(k, ρ, Cₚ)  # Thermal diffusivity [m²/s]
        r = α * Δt / (2 * Δz^2)

        # Initialize the tridiagonal matrix coefficients
        state.solver_cache.a .= -r
        state.solver_cache.a[begin] = 0
        state.solver_cache.a[end]   = 0

        state.solver_cache.b .= 1 + 2r
        state.solver_cache.b[begin] = 1
        state.solver_cache.b[end]   = 1

        state.solver_cache.c .= -r
        state.solver_cache.c[begin] = 0
        state.solver_cache.c[end]   = 0

        # Set the right-hand side vector
        for i_depth in 2:n_depth-1
            state.solver_cache.d[i_depth] = r*T[i_depth+1, i_face] + (1-2r)*T[i_depth, i_face] + r*T[i_depth-1, i_face]
        end
        state.solver_cache.d[begin] = T[begin, i_face]
        state.solver_cache.d[end]   = T[end, i_face]
            
        # Apply lower boundary condition to the matrix system
        if state.problem.lower_boundary_condition isa IsothermalBoundaryCondition
            # T[end] = T_iso
            state.solver_cache.a[end] = 0
            state.solver_cache.b[end] = 1
            state.solver_cache.c[end] = 0
            state.solver_cache.d[end] = state.problem.lower_boundary_condition.T_iso
        elseif state.problem.lower_boundary_condition isa InsulationBoundaryCondition
            # ∂T/∂z = 0 → T[n+1] = T[n]
            state.solver_cache.a[end] = -r
            state.solver_cache.b[end] = 1 + r
            state.solver_cache.c[end] = 0
            # For Crank-Nicolson, modify RHS
            state.solver_cache.d[end] = r*T[end-1, i_face] + (1-r)*T[end, i_face]
        end
            
        # Apply upper boundary condition to the matrix system
        if state.problem.upper_boundary_condition isa IsothermalBoundaryCondition
            # T[1] = T_iso
            state.solver_cache.a[begin] = 0
            state.solver_cache.b[begin] = 1
            state.solver_cache.c[begin] = 0
            state.solver_cache.d[begin] = state.problem.upper_boundary_condition.T_iso
        elseif state.problem.upper_boundary_condition isa InsulationBoundaryCondition
            # ∂T/∂z = 0 → T[0] = T[1]
            state.solver_cache.a[begin] = 0
            state.solver_cache.b[begin] = 1 + r
            state.solver_cache.c[begin] = -r
            # For Crank-Nicolson, modify RHS
            state.solver_cache.d[begin] = r*T[begin+1, i_face] + (1-r)*T[begin, i_face]
        elseif state.problem.upper_boundary_condition isa RadiationBoundaryCondition
            # For radiation BC, we keep the standard interior equation
            # and handle it separately after solving
        end
            
        # Solve the tridiagonal system
        tridiagonal_matrix_algorithm!(state)
            
        # Special handling for radiation boundary condition
        if state.problem.upper_boundary_condition isa RadiationBoundaryCondition
            update_upper_temperature!(state, i_face)
        end
            
        # Copy temperature at next time step
        T[:, i_face] .= state.solver_cache.x
    end
end


"""
    tridiagonal_matrix_algorithm!(a, b, c, d, x)
    tridiagonal_matrix_algorithm!(state::SingleAsteroidThermoPhysicalModel)

Tridiagonal matrix algorithm to solve the heat conduction equation
by the implicit (backward) Euler and Crank-Nicolson methods.

    | b₁ c₁ 0  ⋯  0   | | x₁ |   | d₁ |
    | a₂ b₂ c₂ ⋯  0   | | x₂ |   | d₂ |
    | 0  a₃ b₃ ⋯  0   | | x₃ | = | d₃ |
    | ⋮  ⋮  ⋮  ⋱  cₙ₋₁| | ⋮  |   | ⋮  |
    | 0  0  0  aₙ bₙ  | | xₙ |   | dₙ |     

# References
- https://en.wikipedia.org/wiki/Tridiagonal_matrix_algorithm
"""
function tridiagonal_matrix_algorithm!(a, b, c, d, x)
    N = length(d)

    # Forward sweep
    for i in 2:N
        w = a[i] / b[i - 1]
        b[i] -= w * c[i-1]
        d[i] -= w * d[i-1]
    end

    # Back substitution
    x[N] = d[N] / b[N]
    for i in N-1:-1:1
        x[i] = (d[i] - c[i] * x[i+1]) / b[i]
    end
end

tridiagonal_matrix_algorithm!(state::SingleAsteroidThermoPhysicalState) = tridiagonal_matrix_algorithm!(state.solver_cache.a, state.solver_cache.b, state.solver_cache.c, state.solver_cache.d, state.solver_cache.x)


# ╔═══════════════════════════════════════════════════════════════════╗
# ║                    Upper boundary condition                       ║
# ╚═══════════════════════════════════════════════════════════════════╝

"""
    update_upper_temperature!(state::SingleAsteroidThermoPhysicalState, i::Integer)

Update the temperature of the upper surface based on the boundary condition `state.problem.upper_boundary_condition`.

# Arguments
- `state` : Thermophysical simulation state for a single asteroid
- `i`    : Index of the face of the shape model
"""
function update_upper_temperature!(state::SingleAsteroidThermoPhysicalState, i::Integer)

    #### Radiation boundary condition ####
    if state.problem.upper_boundary_condition isa RadiationBoundaryCondition
        k     = state.problem.thermo_params.thermal_conductivity[i]
        ρ     = state.problem.thermo_params.density[i]
        Cₚ    = state.problem.thermo_params.heat_capacity[i]
        R_vis = state.problem.thermo_params.reflectance_vis[i]
        R_ir  = state.problem.thermo_params.reflectance_ir[i]
        ε     = state.problem.thermo_params.emissivity[i]
        Δz    = state.problem.thermo_params.Δz
    
        F_sun = state.flux_sun[i]
        F_scat = state.flux_scat[i]
        F_rad = state.flux_rad[i]
        F_abs = absorbed_energy_flux(R_vis, R_ir, F_sun, F_scat, F_rad)
        update_surface_temperature!(state.solver_cache.x, F_abs, k, ρ, Cₚ, ε, Δz)
    #### Insulation boundary condition ####
    elseif state.problem.upper_boundary_condition isa InsulationBoundaryCondition
        state.solver_cache.x[begin] = state.solver_cache.x[begin+1]
    #### Isothermal boundary condition ####
    elseif state.problem.upper_boundary_condition isa IsothermalBoundaryCondition
        state.solver_cache.x[begin] = state.problem.upper_boundary_condition.T_iso
    else
        error("The given upper boundary condition is not implemented.")
    end
end


"""
    update_surface_temperature!(T::AbstractVector, F_abs::Real, k::Real, ρ::Real, Cₚ::Real, ε::Real, Δz::Real)

Newton's method to update the surface temperature under radiation boundary condition.

# Arguments
- `T`       : 1-D array of temperatures
- `F_abs`   : Total energy flux absorbed by the facet
- `k`       : Thermal conductivity [W/m/K]
- `ρ`       : Density [kg/m³]
- `Cₚ`      : Heat capacity [J/kg/K]
- `ε`       : Emissivity [-]
- `Δz`      : Depth step width [m]
"""
function update_surface_temperature!(T::AbstractVector, F_abs::Float64, k::Float64, ρ::Float64, Cₚ::Float64, ε::Float64, Δz::Float64)
    εσ = ε * σ_SB  # Pre-compute emissivity × Stefan-Boltzmann constant

    # Newton-Raphson iteration to solve the nonlinear energy balance equation
    for _ in 1:20
        T_pri = T[begin]  # Store previous temperature for convergence check

        # Energy balance at surface: F_absorbed + F_conduction - F_emission = 0
        # f(T) = F_abs + k(T₂-T₁)/Δz - εσT₁⁴ = 0
        f = F_abs + k * (T[begin+1] - T[begin]) / Δz - εσ*T[begin]^4
        
        # Derivative: df/dT = -k/Δz - 4εσT³
        df = - k / Δz - 4*εσ*T[begin]^3             
        
        # Newton-Raphson update: T_new = T_old - f/df
        T[begin] -= f / df

        # Check relative convergence
        err = abs(1 - T_pri / T[begin])
        err < 1e-10 && return
    end
end


# ╔═══════════════════════════════════════════════════════════════════╗
# ║                    Lower boundary condition                       ║
# ╚═══════════════════════════════════════════════════════════════════╝

"""
    update_lower_temperature!(state::SingleAsteroidThermoPhysicalState)

Update the temperature at the lower boundary (deepest layer) based on the boundary condition.

# Arguments
- `state::SingleAsteroidThermoPhysicalState` : Thermophysical simulation state for a single asteroid

# Boundary Conditions
The function applies one of the following boundary conditions at the bottom of the computational domain:

1. **Insulation (Neumann)**: `∂T/∂z = 0`
   - No heat flux through the lower boundary
   - Temperature gradient is zero: `T[end] = T[end-1]`
   - Most commonly used for asteroid modeling

2. **Isothermal (Dirichlet)**: `T = T_iso`
   - Fixed temperature at the lower boundary
   - Used when deep interior temperature is known
   - `T[end] = state.problem.lower_boundary_condition.T_iso`

# Notes
- This function is called after solving the heat conduction equation
- For explicit Euler method, it directly updates the temperature vector
- The lower boundary should be deep enough that the chosen condition doesn't affect surface temperatures
"""
function update_lower_temperature!(state::SingleAsteroidThermoPhysicalState)

    #### Insulation boundary condition ####
    if state.problem.lower_boundary_condition isa InsulationBoundaryCondition
        state.solver_cache.x[end] = state.solver_cache.x[end-1]
    #### Isothermal boundary condition ####
    elseif state.problem.lower_boundary_condition isa IsothermalBoundaryCondition
        state.solver_cache.x[end] = state.problem.lower_boundary_condition.T_iso
    else
        error("The lower boundary condition is not implemented.")
    end
end
