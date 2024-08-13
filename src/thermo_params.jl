

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
- `Γ` : Thermal inertia [tiu (thermal intertia unit)]
"""
thermal_inertia(k, ρ, Cₚ) = @. √(k * ρ * Cₚ)


# ****************************************************************
#               Struct for thermophysical parameters
# ****************************************************************

abstract type AbstractThermoParams end

"""
    struct NonUniformThermoParams

# Fields
- `period`    : Cycle of thermal cycle (rotation period) [sec]
- `skindepth` : Thermal skin depth [m]
- `inertia` : Thermal inertia [tiu (thermal intertia unit)]
- `reflectance_vis` : reflectances for visible light [-]
- `reflectance_ir`  : Reflectances for thermal infrared [-]
- `emissivity` : Emissivity [-]

- `z_max` : Depth of the bottom of a heat conduction equation [m]
- `Δz`    : Depth step width [m]
- `n_depth` : Number of depth steps
"""
struct NonUniformThermoParams <: AbstractThermoParams
    period  ::Float64          # Common for all faces
    skindepth::Vector{Float64}
    inertia ::Vector{Float64}
    reflectance_vis   ::Vector{Float64}
    reflectance_ir    ::Vector{Float64}
    emissivity::Vector{Float64}

    z_max   ::Float64          # Common for all faces
    Δz      ::Float64          # Common for all faces
    n_depth ::Int              # Common for all faces
end

"""
    struct UniformThermoParams

# Fields
- `period`: Thermal cycle (rotation period) [sec]
- `skindepth`: Thermal skin depth [m]
- `inertia`  : Thermal inertia [tiu (thermal intertia unit)]
- `reflectance_vis` : Reflectance for visible light [-]
- `reflectance_ir`  : reflectances for thermal infrared [-]
- `emissivity` : Emissivity [-]

- `z_max` : Depth of the bottom of a heat conduction equation [m]
- `Δz`    : Depth step width [m]
- `n_depth`: Number of depth steps
"""
struct UniformThermoParams <: AbstractThermoParams
    period  ::Float64
    skindepth::Float64
    inertia ::Float64
    reflectance_vis ::Float64
    reflectance_ir  ::Float64
    emissivity::Float64

    z_max   ::Float64
    Δz      ::Float64
    n_depth ::Int
end


"""
    thermoparams(; P, l, Γ, R_vis, R_ir, ε, z_max, n_depth)
"""
function thermoparams(; P, l, Γ, R_vis, R_ir, ε, z_max, n_depth)

    Δz = z_max / (n_depth - 1)
    LENGTH = maximum(length, [R_vis, R_ir, ε, l, Γ])

    if LENGTH > 1
        R_vis   isa Real && (R_vis  = fill(R_vis,  LENGTH))
        R_ir  isa Real && (R_ir = fill(R_ir, LENGTH))
        ε     isa Real && (ε    = fill(ε,    LENGTH))
        l     isa Real && (l    = fill(l,    LENGTH))
        Γ     isa Real && (Γ    = fill(Γ,    LENGTH))
        
        NonUniformThermoParams(P, l, Γ, R_vis, R_ir, ε, z_max, Δz, n_depth)
    else
        UniformThermoParams(P, l, Γ, R_vis, R_ir, ε, z_max, Δz, n_depth)
    end
end


function Base.show(io::IO, params::UniformThermoParams)

    msg =  "⋅-----------------------------------⋅\n"
    msg *= "|     Thermophysical parameters     |\n"
    msg *= "⋅-----------------------------------⋅\n"

    msg *= "  P       = $(params.period) [sec]\n"
    msg *= "          = $(SPICE.convrt(params.period, "seconds", "hours")) [h]\n"
    msg *= "  l       = $(params.skindepth) [m]\n"
    msg *= "  Γ       = $(params.inertia) [tiu]\n"
    msg *= "  R_vis   = $(params.reflectance_vis) [-]\n"
    msg *= "  R_ir    = $(params.reflectance_ir) [-]\n"
    msg *= "  ε       = $(params.emissivity) [-]\n"
  
    msg *= "-----------------------------------\n"

    msg *= "  z_max   = $(params.z_max) [m]\n"
    msg *= "          = $(params.z_max / params.skindepth) [l]\n"
    msg *= "  Δz      = $(params.Δz) [m]\n"
    msg *= "          = $(params.Δz / params.skindepth) [l]\n"
    msg *= "  n_depth = $(params.n_depth)\n"
    
    msg *= "-----------------------------------\n"
    
    print(io, msg)
end

