

# ****************************************************************
#                      1D heat conduction
# ****************************************************************

"""
    update_temperature!(stpm::SingleTPM, Δt)

Calculate the temperature for the next time step based on 1D heat conduction equation.

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
    Nz = size(T, 1)
    Ns = size(T, 2)

    for nₛ in 1:Ns
        P  = stpm.thermo_params.P
        Δz = stpm.thermo_params.Δz
        l  = (stpm.thermo_params.l isa Real ? stpm.thermo_params.l : stpm.thermo_params.l[nₛ])

        λ = (Δt/P) / (Δz/l)^2 / 4π
        λ ≥ 0.5 && error("The forward Euler method is unstable because λ = $λ. This should be less than 0.5.")

        for nz in 2:(Nz-1)
            stpm.SOLVER.T[nz] = (1-2λ)*T[nz, nₛ] + λ*(T[nz+1, nₛ] + T[nz-1, nₛ])  # Predict temperature at next time step
        end

        ## Apply boundary conditions
        update_upper_temperature!(stpm, nₛ)
        update_lower_temperature!(stpm)

        T[:, nₛ] .= stpm.SOLVER.T  # Copy temperature at next time step
    end
end


"""
    backward_euler!(stpm::SingleTPM, Δt)

Predict the temperature at the next time step by the backward Euler method.
- Implicit in time (Unconditionally stable in the heat conduction equation)
- First order in time
- Second order in space
In this function, the heat conduction equation is non-dimensionalized in time and length.
"""
function backward_euler!(stpm::SingleTPM, Δt)
    # T = stpm.temperature
    # Nz = size(T, 1)
    # Ns = size(T, 2)

    # for nₛ in 1:Ns
    #     λ = (stpm.thermo_params.λ isa Real ? stpm.thermo_params.λ : stpm.thermo_params.λ[nₛ])

    #     stpm.SOLVER.a .= -λ
    #     stpm.SOLVER.a[begin] = 0
    #     stpm.SOLVER.a[end]   = 0

    #     stpm.SOLVER.b .= 1 + 2λ
    #     stpm.SOLVER.b[begin] = 1
    #     stpm.SOLVER.b[end]   = 1

    #     stpm.SOLVER.c .= -λ
    #     stpm.SOLVER.c[begin] = 0
    #     stpm.SOLVER.c[end]   = 0

    #     stpm.SOLVER.d .= T[:, nₛ, nₜ]

    #     tridiagonal_matrix_algorithm!(stpm)
    #     T[:, nₛ, nₜ+1] .= stpm.SOLVER.x
    # end

    ## Apply boundary conditions
    # update_upper_temperature!(stpm)
    # update_lower_temperature!(stpm)
end


"""
    crank_nicolson!(stpm::SingleTPM, Δt)

Predict the temperature at the next time step by the Crank-Nicolson method.
- Implicit in time (Unconditionally stable in the heat conduction equation)
- Second order in time
- Second order in space
In this function, the heat conduction equation is non-dimensionalized in time and length.
"""
function crank_nicolson!(stpm::SingleTPM, Δt)
    # T = stpm.temperature
    # Nz = size(T, 1)
    # Ns = size(T, 2)

    # Δt̄ = stpm.thermo_params.Δt / stpm.thermo_params.P  # Non-dimensional timestep, normalized by period
    # Δz̄ = stpm.thermo_params.Δz / stpm.thermo_params.l  # Non-dimensional step in depth, normalized by thermal skin depth
    # r = (1/4π) * (Δt̄ / 2Δz̄^2)

    # for nₛ in 1:Ns
    #     stpm.SOLVER.a .= -r
    #     stpm.SOLVER.a[begin] = 0
    #     stpm.SOLVER.a[end]   = 0

    #     stpm.SOLVER.b .= 1 + 2r
    #     stpm.SOLVER.b[begin] = 1
    #     stpm.SOLVER.b[end]   = 1

    #     stpm.SOLVER.c .= -r
    #     stpm.SOLVER.c[begin] = 0
    #     stpm.SOLVER.c[end]   = 0

    #     for nz in 2:Nz-1
    #         stpm.SOLVER.d[nz] = r*T[nz+1, nₛ, nₜ] + (1-2r)*T[nz, nₛ, nₜ] + r*T[nz-1, nₛ, nₜ]
    #     end

    #     # stpm.SOLVER.d[1]  = 0  # Upper boundary condition
    #     # stpm.SOLVER.d[Nz] = 0 # Lower boundary condition

    #     tridiagonal_matrix_algorithm!(stpm)
    #     T[:, nₛ, nₜ+1] .= stpm.SOLVER.x
    # end

    # ## Apply boundary conditions
    # update_upper_temperature!(stpm)
    # update_lower_temperature!(stpm)
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
    update_upper_temperature!(stpm::SingleTPM, nₛ::Integer)

Update the temperature of the upper surface based on the boundary condition `stpm.BC_UPPER`.

# Arguments
- `stpm`      : Thermophysical model for a single asteroid
- `nₛ`        : Index of the face of the shape model
"""
function update_upper_temperature!(stpm::SingleTPM, nₛ::Integer)

    #### Radiation boundary condition ####
    if stpm.BC_UPPER isa RadiationBoundaryCondition
        P    = stpm.thermo_params.P
        l    = (stpm.thermo_params.l    isa Real ? stpm.thermo_params.l    : stpm.thermo_params.l[nₛ]   )
        Γ    = (stpm.thermo_params.Γ    isa Real ? stpm.thermo_params.Γ    : stpm.thermo_params.Γ[nₛ]   )
        A_B  = (stpm.thermo_params.A_B  isa Real ? stpm.thermo_params.A_B  : stpm.thermo_params.A_B[nₛ] )
        A_TH = (stpm.thermo_params.A_TH isa Real ? stpm.thermo_params.A_TH : stpm.thermo_params.A_TH[nₛ])
        ε    = (stpm.thermo_params.ε    isa Real ? stpm.thermo_params.ε    : stpm.thermo_params.ε[nₛ]   )
        Δz   = stpm.thermo_params.Δz
    
        F_sun, F_scat, F_rad = stpm.flux[nₛ, :]
        F_total = flux_total(A_B, A_TH, F_sun, F_scat, F_rad)
        update_surface_temperature!(stpm.SOLVER.T, F_total, P, l, Γ, ε, Δz)
    #### Insulation boundary condition ####
    elseif stpm.BC_UPPER isa InsulationBoundaryCondition
        stpm.SOLVER.T[begin] = stpm.SOLVER.T[begin+1]
    #### Isothermal boundary condition ####
    elseif stpm.BC_UPPER isa IsothermalBoundaryCondition
        stpm.SOLVER.T[begin] = stpm.BC_UPPER.T_iso
    else
        error("The upper boundary condition is not implemented.")
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
- `ε`       : Emissivity
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
        stpm.SOLVER.T[end] = stpm.SOLVER.T[end-1]
    #### Isothermal boundary condition ####
    elseif stpm.BC_LOWER isa IsothermalBoundaryCondition
        stpm.SOLVER.T[end] = stpm.BC_LOWER.T_iso
    else
        error("The lower boundary condition is not implemented.")
    end
end

