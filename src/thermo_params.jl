

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
- `period`          : Period of thermal cycle (rotation period) [sec]
- `skindepth`       : Vector of thermal skin depth for each facet [m]
- `inertia`         : Vector of thermal inertia for each facet [thermal inertia unit (tiu)]
- `reflectance_vis` : Vector of reflectance in visible light for each facet [-]
- `reflectance_ir`  : Vector of reflectance in thermal infrared for each facet [-]
- `emissivity`      : Vector of emissivity for each facet [-]

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

