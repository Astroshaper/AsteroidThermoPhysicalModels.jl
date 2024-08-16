

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
    struct ThermoParams

# Fields
- `period`          : Cycle of thermal cycle (rotation period) [sec]
- `skindepth`       : Thermal skin depth for each facet [m]
- `inertia`         : Thermal inertia for each facet [thermal inertia unit (tiu)]
- `reflectance_vis` : Reflectance in visible light for each facet [-]
- `reflectance_ir`  : Reflectance in thermal infrared for each facet [-]
- `emissivity`      : Emissivity for each facet [-]

- `z_max`   : Depth of the lower boundary of a heat conduction equation [m]
- `Δz`      : Depth step width [m]
- `n_depth` : Number of depth steps
"""
struct ThermoParams <: AbstractThermoParams
    period          ::Float64
    skindepth       ::Vector{Float64}
    inertia         ::Vector{Float64}
    reflectance_vis ::Vector{Float64}
    reflectance_ir  ::Vector{Float64}
    emissivity      ::Vector{Float64}

    z_max   ::Float64
    Δz      ::Float64
    n_depth ::Int
end


# function thermoparams(;
#     P,
#     l::Vector{Float64},
#     Γ::Vector{Float64},
#     R_vis::Vector{Float64},
#     R_ir::Vector{Float64},
#     ε::Vector{Float64},
#     z_max,
#     n_depth
# )
#     Δz = z_max / (n_depth - 1)    
#     return ThermoParams(P, l, Γ, R_vis, R_ir, ε, z_max, Δz, n_depth)
# end


# function thermoparams(;
#     P,
#     l::Float64,
#     Γ::Float64,
#     R_vis::Float64,
#     R_ir::Float64,
#     ε::Float64,
#     z_max,
#     n_depth
# )
#     Δz = z_max / (n_depth - 1)    
#     return ThermoParams(P, [l], [Γ], [R_vis], [R_ir], [ε], z_max, Δz, n_depth)
# end


# function Base.show(io::IO, params::ThermoParams)

#     msg =  "⋅-----------------------------------⋅\n"
#     msg *= "|     Thermophysical parameters     |\n"
#     msg *= "⋅-----------------------------------⋅\n"

#     msg *= "  Rotation period        = $(params.period) [sec]\n"
#     msg *= "                         = $(SPICE.convrt(params.period, "seconds", "hours")) [h]\n"
#     msg *= "  Thermal skindepth      = $(params.skindepth) [m]\n"
#     msg *= "  Thermal inertia        = $(params.inertia) [tiu]\n"
#     msg *= "  Reflectance (visible)  = $(params.reflectance_vis)\n"
#     msg *= "  Reflectance (infrared) = $(params.reflectance_ir)\n"
#     msg *= "  Emissivity             = $(params.emissivity)\n"
  
#     msg *= "-----------------------------------\n"

#     msg *= "  Maximum depth       = $(params.z_max) [m]\n"
#     msg *= "                      = $(params.z_max / params.skindepth) [l]\n"
#     msg *= "  Depth step width    = $(params.Δz) [m]\n"
#     msg *= "                      = $(params.Δz / params.skindepth) [l]\n"
#     msg *= "  Num. of depth nodes = $(params.n_depth)\n"
    
#     msg *= "-----------------------------------\n"
    
#     print(io, msg)
# end

