

# ****************************************************************
#              Thermal skin depth & Thermal inertia
# ****************************************************************

"""
    thermal_skin_depth(params)      -> l_2π
    thermal_skin_depth(P, k, ρ, Cp) -> l_2π

# Arguments
- `P`  : Cycle of thermal cycle [sec]
- `k`  : Thermal conductivity [W/m/K]
- `ρ`  : Material density [kg/m³]
- `Cₚ` : Heat capacity [J/kg/K]

# Return
- `l_2π` : Thermal skin depth [m], as defined in Rozitis & Green (2011).
"""
thermal_skin_depth(P, k, ρ, Cₚ) = @. √(4π * P * k / (ρ * Cₚ))


"""
    thermal_inertia(params)   -> Γ
    thermal_inertia(k, ρ, Cp) -> Γ

# Arguments
- `k`  : Thermal conductivity [W/m/K]
- `ρ`  : Material density [kg/m³]
- `Cₚ` : Heat capacity [J/kg/K]

# Return
- `Γ` : Thermal inertia [J ⋅ m⁻² ⋅ K⁻¹ ⋅ s⁻⁰⁵ (tiu)]
"""
thermal_inertia(k, ρ, Cₚ) = @. √(k * ρ * Cₚ)


# ****************************************************************
#               Struct for thermophysical properties
# ****************************************************************

abstract type AbstractThermoParams end

"""
    struct NonUniformThermoParams

# Fields
- `P`     : Cycle of thermal cycle (rotation period) [sec]
- `l`     : Thermal skin depth [m]
- `Γ`     : Thermal inertia [J ⋅ m⁻² ⋅ K⁻¹ ⋅ s⁻⁰⁵ (tiu)]
- `λ`     : Non-dimensional coefficient for heat diffusion equation
- `A_B`   : Bond albedo
- `A_TH`  : Albedo at thermal radiation wavelength
- `ε`     : Emissivity

- `t_begin` : Start time of the thermophysical simulation [sec]
- `t_end`   : End time of the thermophysical simulation [sec]
- `Δt`      : Timestep [sec]
- `Nt`      : Number of timesteps

- `z_max` : Depth of the bottom of a heat conduction equation [m]
- `Δz`    : Depth step width [m]
- `Nz`    : Number of depth steps
"""
struct NonUniformThermoParams <: AbstractThermoParams
    P       ::Float64          # Common for all faces
    l       ::Vector{Float64}
    Γ       ::Vector{Float64}
    λ       ::Vector{Float64}
    A_B     ::Vector{Float64}
    A_TH    ::Vector{Float64}
    ε       ::Vector{Float64}

    t_begin ::Float64          # Common for all faces
    t_end   ::Float64          # Common for all faces
    Δt      ::Float64          # Common for all faces
    Nt      ::Int              # Common for all faces

    z_max   ::Float64          # Common for all faces
    Δz      ::Float64          # Common for all faces
    Nz      ::Int              # Common for all faces
end

"""
    struct UniformThermoParams

# Fields
- `P`     : Thermal cycle (rotation period) [sec]
- `l`     : Thermal skin depth [m]
- `Γ`     : Thermal inertia [J ⋅ m⁻² ⋅ K⁻¹ ⋅ s⁻⁰⁵ (tiu)]
- `λ`     : Non-dimensional coefficient for heat diffusion equation
- `A_B`   : Bond albedo
- `A_TH`  : Albedo at thermal radiation wavelength
- `ε`     : Emissivity

- `t_begin` : Start time of the thermophysical simulation [sec]
- `t_end`   : End time of the thermophysical simulation [sec]
- `Δt`      : Timestep [sec]
- `Nt`      : Number of timesteps

- `z_max` : Depth of the bottom of a heat conduction equation [m]
- `Δz`    : Depth step width [m]
- `Nz`    : Number of depth steps
"""
struct UniformThermoParams <: AbstractThermoParams
    P       ::Float64
    l       ::Float64
    Γ       ::Float64
    λ       ::Float64
    A_B     ::Float64
    A_TH    ::Float64
    ε       ::Float64

    t_begin ::Float64
    t_end   ::Float64
    Δt      ::Float64
    Nt      ::Int

    z_max   ::Float64
    Δz      ::Float64
    Nz      ::Int
end


"""
    thermoparams(; A_B, A_TH, k, ρ, Cp, ε, t_begin, t_end, Nt, z_max, Nz, P)
"""
function thermoparams(; P, l, Γ, A_B, A_TH, ε, t_begin, t_end, Nt, z_max, Nz)

    Δt = (t_end - t_begin) / (Nt - 1)
    Δz = z_max / (Nz - 1)

    λ = @. (Δt/P) / (Δz/l)^2 / 4π
    maximum(λ) > 0.5 && error("λ should be smaller than 0.5 for convergence of the forward Euler method.")

    LENGTH = maximum(length, [A_B, A_TH, ε, l, Γ, λ])

    if LENGTH > 1
        A_B   isa Real && (A_B  = fill(A_B,  LENGTH))
        A_TH  isa Real && (A_TH = fill(A_TH, LENGTH))
        ε     isa Real && (ε    = fill(ε,    LENGTH))
        l     isa Real && (l    = fill(l,    LENGTH))
        Γ     isa Real && (Γ    = fill(Γ,    LENGTH))
        λ     isa Real && (λ    = fill(λ,    LENGTH))
        
        NonUniformThermoParams(P, l, Γ, λ, A_B, A_TH, ε, t_begin, t_end, Δt, Nt, z_max, Δz, Nz)
    else
        UniformThermoParams(P, l, Γ, λ, A_B, A_TH, ε, t_begin, t_end, Δt, Nt, z_max, Δz, Nz)
    end
end


function Base.show(io::IO, params::UniformThermoParams)

    msg =  "⋅-----------------------------------⋅\n"
    msg *= "|     Thermophysical parameters     |\n"
    msg *= "⋅-----------------------------------⋅\n"

    msg *= "  P       = $(params.P) [sec]\n"
    msg *= "          = $(SPICE.convrt(params.P, "seconds", "hours")) [h]\n"
    msg *= "  l       = $(params.l) [m]\n"
    msg *= "  Γ       = $(params.Γ) [tiu]\n"
    msg *= "  λ       = $(params.λ)\n"
    msg *= "  A_B     = $(params.A_B)\n"
    msg *= "  A_TH    = $(params.A_TH)\n"
    msg *= "  ε       = $(params.ε)\n"
    
    msg *= "-----------------------------------\n"

    msg *= "  t_begin = $(params.t_begin) [sec]\n"
    msg *= "          = $(params.t_begin / params.P) [P]\n"
    msg *= "  t_end   = $(params.t_end) [sec]\n"
    msg *= "          = $(params.t_end / params.P) [P]\n"
    msg *= "  Δt      = $(params.Δt) [sec]\n"
    msg *= "          = $(params.Δt / params.P) [P]\n"
    msg *= "  Nt      = $(params.Nt)\n"

    msg *= "-----------------------------------\n"

    msg *= "  z_max   = $(params.z_max) [m]\n"
    msg *= "          = $(params.z_max / params.l) [l]\n"
    msg *= "  Δz      = $(params.Δz) [m]\n"
    msg *= "          = $(params.Δz / params.l) [l]\n"
    msg *= "  Nz      = $(params.Nz)\n"
    
    msg *= "-----------------------------------\n"
    
    print(io, msg)
end


# ****************************************************************
#                      1D heat conduction
# ****************************************************************

"""
    forward_temperature(shape::ShapeModel, λ, nₜ::Integer)

Calculate the temperature for the next time step (`nₜ + 1`) based on 1D heat conductivity equation.

TO DO: Allow selection of boundary conditions and solvers

# Arguments
- `shape`  : Shape model
- `params` : Thermophysical parameters
- `nₜ`     : Index of the current time step
"""
function update_temperature!(shape::ShapeModel, params::AbstractThermoParams, nₜ::Integer)
    λ = params.λ
    Tⱼ   = @views shape.temperature[:, :, nₜ  ]
    Tⱼ₊₁ = @views shape.temperature[:, :, nₜ+1]

    ## Forward Euler method
    @. Tⱼ₊₁[begin+1:end-1, :] = @views (1-2λ')*Tⱼ[begin+1:end-1, :] + λ'*(Tⱼ[begin+2:end, :] + Tⱼ[begin:end-2, :])

    ## Boundary conditions
    update_surface_temperature!(shape, params, nₜ+1)  # Radiation (surface)
    update_bottom_temperature!(shape, nₜ+1)           # Insulation (bottom)
end


# ****************************************************************
#                   Surface boundary condition
# ****************************************************************

"""
    update_surface_temperature!(shape::ShapeModel, params::AbstractThermoParams, nₜ::Integer)

Update surface temperature under radiative boundary condition using Newton's method

# Arguments
- `shape`  : Shape model (`ShapeModel`)
- `params` : Thermophysical prameters
- `nₜ`     : Index of the current time step
"""
function update_surface_temperature!(shape::ShapeModel, params::AbstractThermoParams, nₜ::Integer)
    for i in eachindex(shape.faces)
        P    = params.P
        l    = (params.l    isa Real ? params.l    : params.l[i]   )
        Γ    = (params.Γ    isa Real ? params.Γ    : params.Γ[i]   )
        A_B  = (params.A_B  isa Real ? params.A_B  : params.A_B[i] )
        A_TH = (params.A_TH isa Real ? params.A_TH : params.A_TH[i])
        ε    = (params.ε    isa Real ? params.ε    : params.ε[i]   )
        Δz   = params.Δz

        F_sun, F_scat, F_rad = shape.flux[i, :]
        F_total = total_flux(A_B, A_TH, F_sun, F_scat, F_rad)
        update_surface_temperature!((@views shape.temperature[:, i, nₜ]), F_total, P, l, Γ, ε, Δz)  # Δz should be normalized by l
    end
end


"""
    update_surface_temperature!(T::AbstractVector, F_total::Real, k::Real, l::Real, Δz::Real, ε::Real)

Newton's method to update the surface temperature under radiative boundary condition.

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
#                   Bottom boundary condition
# ****************************************************************

"""
    update_bottom_temperature!(shape::ShapeModel, nₜ::Integer)

Update bottom temperature under boundary condition of insulation
"""
function update_bottom_temperature!(shape::ShapeModel, nₜ::Integer)
    for i in eachindex(shape.faces)
        shape.temperature[end, i, nₜ] = shape.temperature[end-1, i, nₜ]
    end
end

