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
    update_temperature_zero_conductivity!(stpm::SingleAsteroidTPM)

Update surface temperature for the zero thermal conductivity case.
When thermal conductivity is zero, there is no heat conduction into the subsurface,
and the surface temperature is determined solely by instantaneous radiative equilibrium.

# Arguments
- `stpm::SingleAsteroidTPM` : Thermophysical model for a single asteroid

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
- This function is called when `stpm.thermo_params.thermal_conductivity` is zero
- The temperature instantly adjusts to balance incoming and outgoing radiation
- No subsurface temperatures are updated (only surface layer)
"""
function update_temperature_zero_conductivity!(stpm::SingleAsteroidTPM)
    for i_face in eachindex(stpm.shape.faces)
        R_vis = stpm.thermo_params.reflectance_vis[i_face]
        R_ir  = stpm.thermo_params.reflectance_ir[i_face]
        ε     = stpm.thermo_params.emissivity[i_face]
        εσ    = ε * σ_SB

        F_sun = stpm.flux_sun[i_face]
        F_scat = stpm.flux_scat[i_face]
        F_rad = stpm.flux_rad[i_face]
        F_abs = absorbed_energy_flux(R_vis, R_ir, F_sun, F_scat, F_rad)

        stpm.temperature[begin, i_face] = (F_abs / εσ)^(1/4)
    end
end

"""
    update_temperature!(stpm::SingleAsteroidTPM, Δt)

Update the temperature distribution for the next time step by solving the 1D heat conduction equation.
The solver method is determined by `stpm.SOLVER`, and special handling is applied for zero conductivity.

# Arguments
- `stpm::SingleAsteroidTPM` : Thermophysical model for a single asteroid
- `Δt::Real` : Time step [s]

# Solver Selection
The function automatically selects the appropriate solver based on `stpm.SOLVER`:
- `ExplicitEulerSolver`: Forward Euler method (conditionally stable, requires λ < 0.5)
- `ImplicitEulerSolver`: Backward Euler method (unconditionally stable)
- `CrankNicolsonSolver`: Crank-Nicolson method (unconditionally stable, second-order accurate)

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
function update_temperature!(stpm::SingleAsteroidTPM, Δt)
    # Handle zero-conductivity case
    if iszero(stpm.thermo_params.thermal_conductivity)
        update_temperature_zero_conductivity!(stpm)
        return
    end
    
    # Non-zero conductivity: use selected solver
    if stpm.SOLVER isa ExplicitEulerSolver
        explicit_euler!(stpm, Δt)
    elseif stpm.SOLVER isa ImplicitEulerSolver
        implicit_euler!(stpm, Δt)
    elseif stpm.SOLVER isa CrankNicolsonSolver
        crank_nicolson!(stpm, Δt)
    else
        error("Unknown solver type: $(typeof(stpm.SOLVER)). Expected ExplicitEulerSolver, ImplicitEulerSolver, or CrankNicolsonSolver.")
    end
end


"""
    update_temperature!(btpm::BinaryAsteroidTPM, Δt)

Calculate the temperature for the next time step based on 1D heat conductivity equation.

# Arguments
- `btpm` : Thermophysical model for a binary asteroid
- `Δt`   : Time step [sec]
"""
function update_temperature!(btpm::BinaryAsteroidTPM, Δt)
    update_temperature!(btpm.pri, Δt)
    update_temperature!(btpm.sec, Δt)
end

# ╔═══════════════════════════════════════════════════════════════════╗
# ║               Solvers of a heat conduction equation               ║
# ╚═══════════════════════════════════════════════════════════════════╝


"""
    explicit_euler!(stpm::SingleAsteroidTPM, Δt)

Solve the 1D heat conduction equation using the explicit (forward) Euler method.
This method is conditionally stable and requires careful time step selection.

# Arguments
- `stpm::SingleAsteroidTPM` : Thermophysical model for a single asteroid
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
function explicit_euler!(stpm::SingleAsteroidTPM, Δt)
    T = stpm.temperature
    n_depth = size(T, 1)
    n_face = size(T, 2)

    for i_face in 1:n_face
        k  = stpm.thermo_params.thermal_conductivity[i_face]
        ρ  = stpm.thermo_params.density[i_face]
        Cₚ = stpm.thermo_params.heat_capacity[i_face]
        Δz = stpm.thermo_params.Δz

        α = thermal_diffusivity(k, ρ, Cₚ)  # Thermal diffusivity [m²/s]
        λ = α * Δt / Δz^2
        λ ≥ 0.5 && error("Explicit Euler method is unstable because λ = αΔt/Δz² = $λ (must be < 0.5). Consider reducing time step from Δt = $Δt or using an implicit solver.")

        for i_depth in 2:(n_depth-1)
            stpm.SOLVER.x[i_depth] = (1-2λ)*T[i_depth, i_face] + λ*(T[i_depth+1, i_face] + T[i_depth-1, i_face])  # Predict temperature at next time step
        end

        ## Apply boundary conditions
        update_upper_temperature!(stpm, i_face)
        update_lower_temperature!(stpm)

        T[:, i_face] .= stpm.SOLVER.x  # Copy temperature at next time step
    end
end


"""
    implicit_euler!(stpm::SingleAsteroidTPM, Δt)

Solve the 1D heat conduction equation using the implicit (backward) Euler method.
This method is unconditionally stable, allowing for larger time steps than explicit methods.

# Arguments
- `stpm::SingleAsteroidTPM` : Thermophysical model for a single asteroid
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
function implicit_euler!(stpm::SingleAsteroidTPM, Δt)
    T = stpm.temperature
    n_depth = size(T, 1)
    n_face = size(T, 2)

    for i_face in 1:n_face
        k  = stpm.thermo_params.thermal_conductivity[i_face]
        ρ  = stpm.thermo_params.density[i_face]
        Cₚ = stpm.thermo_params.heat_capacity[i_face]
        Δz = stpm.thermo_params.Δz

        α = thermal_diffusivity(k, ρ, Cₚ)  # Thermal diffusivity [m²/s]
        λ = α * Δt / Δz^2

        # Initialize the tridiagonal matrix coefficients
        stpm.SOLVER.a .= -λ
        stpm.SOLVER.a[begin] = 0
        stpm.SOLVER.a[end]   = 0

        stpm.SOLVER.b .= 1 + 2λ
        stpm.SOLVER.b[begin] = 1
        stpm.SOLVER.b[end]   = 1

        stpm.SOLVER.c .= -λ
        stpm.SOLVER.c[begin] = 0
        stpm.SOLVER.c[end]   = 0

        # Set the right-hand side vector
        stpm.SOLVER.d .= T[:, i_face]
            
        # Apply lower boundary condition to the matrix system
        if stpm.BC_LOWER isa IsothermalBoundaryCondition
            # T[end] = T_iso
            stpm.SOLVER.a[end] = 0
            stpm.SOLVER.b[end] = 1
            stpm.SOLVER.c[end] = 0
            stpm.SOLVER.d[end] = stpm.BC_LOWER.T_iso
        elseif stpm.BC_LOWER isa InsulationBoundaryCondition
            # ∂T/∂z = 0 → T[n+1] = T[n]
            stpm.SOLVER.a[end] = -λ
            stpm.SOLVER.b[end] = 1 + λ
            stpm.SOLVER.c[end] = 0
            # d[end] already contains T[end, i_face]
        end

        # Apply upper boundary condition to the matrix system
        if stpm.BC_UPPER isa IsothermalBoundaryCondition
            # T[1] = T_iso
            stpm.SOLVER.a[begin] = 0
            stpm.SOLVER.b[begin] = 1
            stpm.SOLVER.c[begin] = 0
            stpm.SOLVER.d[begin] = stpm.BC_UPPER.T_iso
        elseif stpm.BC_UPPER isa InsulationBoundaryCondition
            # ∂T/∂z = 0 → T[0] = T[1]
            stpm.SOLVER.a[begin] = 0
            stpm.SOLVER.b[begin] = 1 + λ
            stpm.SOLVER.c[begin] = -λ
            # d[begin] already contains T[begin, i_face]
        elseif stpm.BC_UPPER isa RadiationBoundaryCondition
            # For radiation BC, we keep the standard interior equation
            # and handle it separately after solving
        end
            
        # Solve the tridiagonal system
        tridiagonal_matrix_algorithm!(stpm)
            
        # Special handling for radiation boundary condition
        if stpm.BC_UPPER isa RadiationBoundaryCondition
            update_upper_temperature!(stpm, i_face)
        end
            
        # Copy temperature at next time step
        T[:, i_face] .= stpm.SOLVER.x
    end
end


"""
    crank_nicolson!(stpm::SingleAsteroidTPM, Δt)

Solve the 1D heat conduction equation using the Crank-Nicolson method.
This method combines the explicit and implicit Euler methods for improved accuracy.

# Arguments
- `stpm::SingleAsteroidTPM` : Thermophysical model for a single asteroid
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
function crank_nicolson!(stpm::SingleAsteroidTPM, Δt)
    T = stpm.temperature
    n_depth = size(T, 1)
    n_face = size(T, 2)

    for i_face in 1:n_face
        k  = stpm.thermo_params.thermal_conductivity[i_face]
        ρ  = stpm.thermo_params.density[i_face]
        Cₚ = stpm.thermo_params.heat_capacity[i_face]
        Δz = stpm.thermo_params.Δz

        α = thermal_diffusivity(k, ρ, Cₚ)  # Thermal diffusivity [m²/s]
        r = α * Δt / (2 * Δz^2)

        # Initialize the tridiagonal matrix coefficients
        stpm.SOLVER.a .= -r
        stpm.SOLVER.a[begin] = 0
        stpm.SOLVER.a[end]   = 0

        stpm.SOLVER.b .= 1 + 2r
        stpm.SOLVER.b[begin] = 1
        stpm.SOLVER.b[end]   = 1

        stpm.SOLVER.c .= -r
        stpm.SOLVER.c[begin] = 0
        stpm.SOLVER.c[end]   = 0

        # Set the right-hand side vector
        for i_depth in 2:n_depth-1
            stpm.SOLVER.d[i_depth] = r*T[i_depth+1, i_face] + (1-2r)*T[i_depth, i_face] + r*T[i_depth-1, i_face]
        end
        stpm.SOLVER.d[begin] = T[begin, i_face]
        stpm.SOLVER.d[end]   = T[end, i_face]
            
        # Apply lower boundary condition to the matrix system
        if stpm.BC_LOWER isa IsothermalBoundaryCondition
            # T[end] = T_iso
            stpm.SOLVER.a[end] = 0
            stpm.SOLVER.b[end] = 1
            stpm.SOLVER.c[end] = 0
            stpm.SOLVER.d[end] = stpm.BC_LOWER.T_iso
        elseif stpm.BC_LOWER isa InsulationBoundaryCondition
            # ∂T/∂z = 0 → T[n+1] = T[n]
            stpm.SOLVER.a[end] = -r
            stpm.SOLVER.b[end] = 1 + r
            stpm.SOLVER.c[end] = 0
            # For Crank-Nicolson, modify RHS
            stpm.SOLVER.d[end] = r*T[end-1, i_face] + (1-r)*T[end, i_face]
        end
            
        # Apply upper boundary condition to the matrix system
        if stpm.BC_UPPER isa IsothermalBoundaryCondition
            # T[1] = T_iso
            stpm.SOLVER.a[begin] = 0
            stpm.SOLVER.b[begin] = 1
            stpm.SOLVER.c[begin] = 0
            stpm.SOLVER.d[begin] = stpm.BC_UPPER.T_iso
        elseif stpm.BC_UPPER isa InsulationBoundaryCondition
            # ∂T/∂z = 0 → T[0] = T[1]
            stpm.SOLVER.a[begin] = 0
            stpm.SOLVER.b[begin] = 1 + r
            stpm.SOLVER.c[begin] = -r
            # For Crank-Nicolson, modify RHS
            stpm.SOLVER.d[begin] = r*T[begin+1, i_face] + (1-r)*T[begin, i_face]
        elseif stpm.BC_UPPER isa RadiationBoundaryCondition
            # For radiation BC, we keep the standard interior equation
            # and handle it separately after solving
        end
            
        # Solve the tridiagonal system
        tridiagonal_matrix_algorithm!(stpm)
            
        # Special handling for radiation boundary condition
        if stpm.BC_UPPER isa RadiationBoundaryCondition
            update_upper_temperature!(stpm, i_face)
        end
            
        # Copy temperature at next time step
        T[:, i_face] .= stpm.SOLVER.x
    end
end


"""
    tridiagonal_matrix_algorithm!(a, b, c, d, x)
    tridiagonal_matrix_algorithm!(stpm::SingleAsteroidThermoPhysicalModel)

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

tridiagonal_matrix_algorithm!(stpm::SingleAsteroidTPM) = tridiagonal_matrix_algorithm!(stpm.SOLVER.a, stpm.SOLVER.b, stpm.SOLVER.c, stpm.SOLVER.d, stpm.SOLVER.x)


# ╔═══════════════════════════════════════════════════════════════════╗
# ║                    Upper boundary condition                       ║
# ╚═══════════════════════════════════════════════════════════════════╝

"""
    update_upper_temperature!(stpm::SingleAsteroidTPM, i::Integer)

Update the temperature of the upper surface based on the boundary condition `stpm.BC_UPPER`.

# Arguments
- `stpm` : Thermophysical model for a single asteroid
- `i`    : Index of the face of the shape model
"""
function update_upper_temperature!(stpm::SingleAsteroidTPM, i::Integer)

    #### Radiation boundary condition ####
    if stpm.BC_UPPER isa RadiationBoundaryCondition
        k     = stpm.thermo_params.thermal_conductivity[i]
        ρ     = stpm.thermo_params.density[i]
        Cₚ    = stpm.thermo_params.heat_capacity[i]
        R_vis = stpm.thermo_params.reflectance_vis[i]
        R_ir  = stpm.thermo_params.reflectance_ir[i]
        ε     = stpm.thermo_params.emissivity[i]
        Δz    = stpm.thermo_params.Δz
    
        F_sun = stpm.flux_sun[i]
        F_scat = stpm.flux_scat[i]
        F_rad = stpm.flux_rad[i]
        F_abs = absorbed_energy_flux(R_vis, R_ir, F_sun, F_scat, F_rad)
        update_surface_temperature!(stpm.SOLVER.x, F_abs, k, ρ, Cₚ, ε, Δz)
    #### Insulation boundary condition ####
    elseif stpm.BC_UPPER isa InsulationBoundaryCondition
        stpm.SOLVER.x[begin] = stpm.SOLVER.x[begin+1]
    #### Isothermal boundary condition ####
    elseif stpm.BC_UPPER isa IsothermalBoundaryCondition
        stpm.SOLVER.x[begin] = stpm.BC_UPPER.T_iso
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
    update_lower_temperature!(stpm::SingleAsteroidTPM)

Update the temperature at the lower boundary (deepest layer) based on the boundary condition.

# Arguments
- `stpm::SingleAsteroidTPM` : Thermophysical model for a single asteroid

# Boundary Conditions
The function applies one of the following boundary conditions at the bottom of the computational domain:

1. **Insulation (Neumann)**: `∂T/∂z = 0`
   - No heat flux through the lower boundary
   - Temperature gradient is zero: `T[end] = T[end-1]`
   - Most commonly used for asteroid modeling

2. **Isothermal (Dirichlet)**: `T = T_iso`
   - Fixed temperature at the lower boundary
   - Used when deep interior temperature is known
   - `T[end] = stpm.BC_LOWER.T_iso`

# Notes
- This function is called after solving the heat conduction equation
- For explicit Euler method, it directly updates the temperature vector
- The lower boundary should be deep enough that the chosen condition doesn't affect surface temperatures
"""
function update_lower_temperature!(stpm::SingleAsteroidTPM)

    #### Insulation boundary condition ####
    if stpm.BC_LOWER isa InsulationBoundaryCondition
        stpm.SOLVER.x[end] = stpm.SOLVER.x[end-1]
    #### Isothermal boundary condition ####
    elseif stpm.BC_LOWER isa IsothermalBoundaryCondition
        stpm.SOLVER.x[end] = stpm.BC_LOWER.T_iso
    else
        error("The lower boundary condition is not implemented.")
    end
end
