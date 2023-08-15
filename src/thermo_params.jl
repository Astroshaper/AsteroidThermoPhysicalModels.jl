

# ****************************************************************
#              Thermal skin depth & Thermal inertia
# ****************************************************************

"""
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
#               Struct for thermophysical parameters
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

