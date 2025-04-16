

# ****************************************************************
#                      1D heat conduction
# ****************************************************************

"""
    update_temperature!(stpm::SingleTPM, Δt)

Calculate the temperature for the next time step based on 1D heat conduction equation.
If the thermal inertia (conductivity) is zero, omit to solve the heat conduction equation.
The surface termperature is determined only by radiative equilibrium.

# Arguments
- `stpm` : Thermophysical model for a single asteroid
- `Δt`   : Time step [sec]
"""
function update_temperature!(stpm::SingleTPM, Δt)
    if stpm.SOLVER isa ForwardEulerSolver
        forward_euler!(stpm, Δt)
    elseif stpm.SOLVER isa BackwardEulerSolver
        backward_euler!(stpm, Δt)
    elseif stpm.SOLVER isa CrankNicolsonSolver
        crank_nicolson!(stpm, Δt)
    else
        error("The solver is not implemented.")
    end
end


"""
    update_temperature!(btpm::BinaryTPM, Δt)

Calculate the temperature for the next time step based on 1D heat conductivity equation.

# Arguments
- `btpm` : Thermophysical model for a binary asteroid
- `Δt`   : Time step [sec]
"""
function update_temperature!(btpm::BinaryTPM, Δt)
    update_temperature!(btpm.pri, Δt)
    update_temperature!(btpm.sec, Δt)
end

# ****************************************************************
#               Solvers of a heat conduction equation
# ****************************************************************


"""
    forward_euler!(stpm::SingleTPM, Δt)

Predict the temperature at the next time step by the forward Euler method.
- Explicit in time
- First order in time
In this function, the heat conduction equation is non-dimensionalized in time and length.

# Arguments
- `stpm` : Thermophysical model for a single asteroid
- `Δt`   : Time step [sec]
"""
function forward_euler!(stpm::SingleTPM, Δt)
    T = stpm.temperature
    n_depth = size(T, 1)
    n_face = size(T, 2)

    ## Zero-conductivity (thermal inertia) case
    if iszero(stpm.thermo_params.inertia)
        for i_face in 1:n_face
            R_vis = stpm.thermo_params.reflectance_vis[i_face]
            R_ir  = stpm.thermo_params.reflectance_ir[i_face]
            ε     = stpm.thermo_params.emissivity[i_face]
            εσ = ε * σ_SB

            F_sun, F_scat, F_rad = stpm.flux[i_face, :]
            F_total = flux_total(R_vis, R_ir, F_sun, F_scat, F_rad)

            stpm.temperature[begin, i_face] = (F_total / εσ)^(1/4)
        end
    ## Non-zero-conductivity (thermal inertia) case
    else
        for i_face in 1:n_face
            P  = stpm.thermo_params.period
            Δz = stpm.thermo_params.Δz
            l  = stpm.thermo_params.skindepth[i_face]

            λ = (Δt/P) / (Δz/l)^2 / 4π
            λ ≥ 0.5 && error("The forward Euler method is unstable because λ = $λ. This should be less than 0.5.")

            for i_depth in 2:(n_depth-1)
                stpm.SOLVER.x[i_depth] = (1-2λ)*T[i_depth, i_face] + λ*(T[i_depth+1, i_face] + T[i_depth-1, i_face])  # Predict temperature at next time step
            end

            ## Apply boundary conditions
            update_upper_temperature!(stpm, i_face)
            update_lower_temperature!(stpm)

            T[:, i_face] .= stpm.SOLVER.x  # Copy temperature at next time step
        end
    end
end


"""
    backward_euler!(stpm::SingleTPM, Δt)

Predict the temperature at the next time step by the backward Euler method.
- Implicit in time (Unconditionally stable in the heat conduction equation)
- First order in time
- Second order in space
In this function, the heat conduction equation is non-dimensionalized in time and length.

# Arguments
- `stpm` : Thermophysical model for a single asteroid
- `Δt`   : Time step [sec]
"""
function backward_euler!(stpm::SingleTPM, Δt)
    T = stpm.temperature
    n_depth = size(T, 1)
    n_face = size(T, 2)

    ## Zero-conductivity (thermal inertia) case
    if iszero(stpm.thermo_params.inertia)
        for i_face in 1:n_face
            R_vis = stpm.thermo_params.reflectance_vis[i_face]
            R_ir  = stpm.thermo_params.reflectance_ir[i_face]
            ε     = stpm.thermo_params.emissivity[i_face]
            εσ = ε * σ_SB

            F_sun, F_scat, F_rad = stpm.flux[i_face, :]
            F_total = flux_total(R_vis, R_ir, F_sun, F_scat, F_rad)

            stpm.temperature[begin, i_face] = (F_total / εσ)^(1/4)
        end
    ## Non-zero-conductivity (thermal inertia) case
    else
        for i_face in 1:n_face
            P  = stpm.thermo_params.period
            Δz = stpm.thermo_params.Δz
            l  = stpm.thermo_params.skindepth[i_face]

            # Non-dimensional timestep, normalized by period
            Δt̄ = Δt / P
            # Non-dimensional step in depth, normalized by thermal skin depth
            Δz̄ = Δz / l
            # Stability parameter for implicit methods
            λ = (Δt̄) / (Δz̄^2) / 4π

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
            
            # Apply boundary conditions
            update_lower_temperature!(stpm)
            
            # Solve the tridiagonal system
            tridiagonal_matrix_algorithm!(stpm)
            
            # Apply upper boundary condition (surface temperature)
            stpm.SOLVER.x[begin] = T[begin, i_face]  # 一時的に元の表面温度を保存
            update_upper_temperature!(stpm, i_face)
            
            # Copy temperature at next time step
            T[:, i_face] .= stpm.SOLVER.x
        end
    end
end


"""
    crank_nicolson!(stpm::SingleTPM, Δt)

Predict the temperature at the next time step by the Crank-Nicolson method.
- Implicit in time (Unconditionally stable in the heat conduction equation)
- Second order in time
- Second order in space
In this function, the heat conduction equation is non-dimensionalized in time and length.

# Arguments
- `stpm` : Thermophysical model for a single asteroid
- `Δt`   : Time step [sec]
"""
function crank_nicolson!(stpm::SingleTPM, Δt)
    T = stpm.temperature
    n_depth = size(T, 1)
    n_face = size(T, 2)

    ## Zero-conductivity (thermal inertia) case
    if iszero(stpm.thermo_params.inertia)
        for i_face in 1:n_face
            R_vis = stpm.thermo_params.reflectance_vis[i_face]
            R_ir  = stpm.thermo_params.reflectance_ir[i_face]
            ε     = stpm.thermo_params.emissivity[i_face]
            εσ = ε * σ_SB

            F_sun, F_scat, F_rad = stpm.flux[i_face, :]
            F_total = flux_total(R_vis, R_ir, F_sun, F_scat, F_rad)

            stpm.temperature[begin, i_face] = (F_total / εσ)^(1/4)
        end
    ## Non-zero-conductivity (thermal inertia) case
    else
        for i_face in 1:n_face
            P  = stpm.thermo_params.period
            Δz = stpm.thermo_params.Δz
            l  = stpm.thermo_params.skindepth[i_face]

            # Non-dimensional timestep, normalized by period
            Δt̄ = Δt / P
            # Non-dimensional step in depth, normalized by thermal skin depth
            Δz̄ = Δz / l
            # Stability parameter for Crank-Nicolson method
            r = (1/4π) * (Δt̄ / (2*Δz̄^2))

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
            stpm.SOLVER.d[begin] = T[begin, i_face]  # Will be overwritten by boundary condition
            stpm.SOLVER.d[end]   = T[end, i_face]    # Will be overwritten by boundary condition
            
            # Apply boundary conditions
            update_lower_temperature!(stpm)
            
            # Solve the tridiagonal system
            tridiagonal_matrix_algorithm!(stpm)
            
            # Apply upper boundary condition (surface temperature)
            stpm.SOLVER.x[begin] = T[begin, i_face]  # 一時的に元の表面温度を保存
            update_upper_temperature!(stpm, i_face)
            
            # Copy temperature at next time step
            T[:, i_face] .= stpm.SOLVER.x
        end
    end
end


"""
    tridiagonal_matrix_algorithm!(a, b, c, d, x)
    tridiagonal_matrix_algorithm!(stpm::SingleTPM)

Tridiagonal matrix algorithm to solve the heat conduction equation by the backward Euler and Crank-Nicolson methods.

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

tridiagonal_matrix_algorithm!(stpm::SingleTPM) = tridiagonal_matrix_algorithm!(stpm.SOLVER.a, stpm.SOLVER.b, stpm.SOLVER.c, stpm.SOLVER.d, stpm.SOLVER.x)


# ****************************************************************
#                    Upper boundary condition
# ****************************************************************

"""
    update_upper_temperature!(stpm::SingleTPM, i::Integer)

Update the temperature of the upper surface based on the boundary condition `stpm.BC_UPPER`.

# Arguments
- `stpm` : Thermophysical model for a single asteroid
- `i`    : Index of the face of the shape model
"""
function update_upper_temperature!(stpm::SingleTPM, i::Integer)

    #### Radiation boundary condition ####
    if stpm.BC_UPPER isa RadiationBoundaryCondition
        P     = stpm.thermo_params.period
        l     = stpm.thermo_params.skindepth[i]
        Γ     = stpm.thermo_params.inertia[i]
        R_vis = stpm.thermo_params.reflectance_vis[i]
        R_ir  = stpm.thermo_params.reflectance_ir[i]
        ε     = stpm.thermo_params.emissivity[i]
        Δz    = stpm.thermo_params.Δz
    
        F_sun, F_scat, F_rad = stpm.flux[i, :]
        F_total = flux_total(R_vis, R_ir, F_sun, F_scat, F_rad)
        update_surface_temperature!(stpm.SOLVER.x, F_total, P, l, Γ, ε, Δz)
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
    update_surface_temperature!(T::AbstractVector, F_total::Real, k::Real, l::Real, Δz::Real, ε::Real)

Newton's method to update the surface temperature under radiation boundary condition.

# Arguments
- `T`       : 1-D array of temperatures
- `F_total` : Total energy absorbed by the facet
- `Γ`       : Thermal inertia [tiu]
- `P`       : Period of thermal cycle [sec]
- `Δz̄`      : Non-dimensional step in depth, normalized by thermal skin depth `l`
- `ε`       : Emissivity [-]
"""
function update_surface_temperature!(T::AbstractVector, F_total::Float64, P::Float64, l::Float64, Γ::Float64, ε::Float64, Δz::Float64)
    Δz̄ = Δz / l    # Dimensionless length of depth step
    εσ = ε * σ_SB

    for _ in 1:20
        T_pri = T[begin]

        f = F_total + Γ / √(4π * P) * (T[begin+1] - T[begin]) / Δz̄ - εσ*T[begin]^4
        df = - Γ / √(4π * P) / Δz̄ - 4*εσ*T[begin]^3             
        T[begin] -= f / df

        err = abs(1 - T_pri / T[begin])
        err < 1e-10 && return
    end
end


# ****************************************************************
#                    Lower boundary condition
# ****************************************************************

"""
    update_bottom_temperature!(shape::ShapeModel)

Update the temperature of the bottom surface based on the boundary condition `stpm.BC_LOWER`.

# Arguments
- `stpm`       : Thermophysical model for a single asteroid
"""
function update_lower_temperature!(stpm::SingleTPM)

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
