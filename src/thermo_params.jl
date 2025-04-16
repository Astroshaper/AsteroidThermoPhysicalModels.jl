

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
- `thermal_conductivity` : Vector of thermal conductivity for each facet [W/m/K]
- `density`              : Vector of density for each facet [kg/m³]
- `heat_capacity`        : Vector of heat capacity for each facet [J/kg/K]

- `reflectance_vis` : Vector of reflectance in visible light for each facet [-]
- `reflectance_ir`  : Vector of reflectance in thermal infrared for each facet [-]
- `emissivity`      : Vector of emissivity for each facet [-]

- `z_max`   : Depth of the lower boundary of a heat conduction equation [m]
- `Δz`      : Depth step width [m]
- `n_depth` : Number of depth steps
"""
struct ThermoParams <: AbstractThermoParams
    thermal_conductivity ::Vector{Float64}
    density              ::Vector{Float64}
    heat_capacity        ::Vector{Float64}

    reflectance_vis ::Vector{Float64}
    reflectance_ir  ::Vector{Float64}
    emissivity      ::Vector{Float64}

    z_max   ::Float64
    Δz      ::Float64
    n_depth ::Int
end


"""
    ThermoParams(
        thermal_conductivity ::Float64,
        density              ::Float64,
        heat_capacity        ::Float64,
        reflectance_vis      ::Float64,
        reflectance_ir       ::Float64,
        emissivity           ::Float64,
        z_max                ::Float64,
        Δz                   ::Float64,
        n_depth              ::Int
    )

Outer constructor for `ThermoParams`.
You can give the same parameters to all facets by `Float64`.

# Arguments
- `thermal_conductivity` : Thermal conductivity [W/m/K]
- `density`              : Density [kg/m³]
- `heat_capacity`        : Heat capacity [J/kg/K]

- `reflectance_vis` : Reflectance in visible light [-]
- `reflectance_ir`  : Reflectance in thermal infrared [-]
- `emissivity`      : Emissivity [-]

- `z_max`           : Depth of the lower boundary of a heat conduction equation [m]
- `Δz`              : Depth step width [m]
- `n_depth`         : Number of depth steps
"""
function ThermoParams(
    thermal_conductivity ::Float64,
    density              ::Float64,
    heat_capacity        ::Float64,
    reflectance_vis      ::Float64,
    reflectance_ir       ::Float64,
    emissivity           ::Float64,
    z_max                ::Float64,
    Δz                   ::Float64,
    n_depth              ::Int
)

    return ThermoParams([thermal_conductivity], [density], [heat_capacity], [reflectance_vis], [reflectance_ir], [emissivity], z_max, Δz, n_depth)
end
