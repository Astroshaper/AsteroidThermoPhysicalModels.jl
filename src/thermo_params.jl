

# ****************************************************************
#              Thermal skin depth & Thermal inertia
# ****************************************************************

"""
    thermal_skin_depth(P, k, ρ, Cp) -> l_2π

Calculate the thermal skin depth for a periodic temperature variation.

# Arguments
- `P::Real`  : Period of thermal cycle [s]
- `k::Real`  : Thermal conductivity [W/m/K]
- `ρ::Real`  : Material density [kg/m³]
- `Cₚ::Real` : Heat capacity [J/kg/K]

# Returns
- `l_2π::Real` : Thermal skin depth [m]

# Mathematical Formula
The thermal skin depth is defined as:
```math
l_{2\\pi} = \\sqrt{\\frac{4\\pi P k}{\\rho C_p}}
```

# Physical Meaning
- Represents the e-folding depth of temperature variations
- Temperature amplitude decreases by factor e^(-2π) ≈ 0.0019 at this depth
- Useful for determining computational domain depth

# Reference
- Rozitis & Green (2011), MNRAS 415, 2042-2062
"""
thermal_skin_depth(P, k, ρ, Cₚ) = @. √(4π * P * k / (ρ * Cₚ))


"""
    thermal_inertia(k, ρ, Cp) -> Γ

Calculate the thermal inertia of a material.

# Arguments
- `k::Real`  : Thermal conductivity [W/m/K]
- `ρ::Real`  : Material density [kg/m³]
- `Cₚ::Real` : Heat capacity [J/kg/K]

# Returns
- `Γ::Real` : Thermal inertia [J m⁻² K⁻¹ s⁻¹/²]

# Mathematical Formula
```math
\\Gamma = \\sqrt{k \\rho C_p}
```

# Physical Meaning
- Measures resistance to temperature change
- High Γ: slow temperature response (rock-like)
- Low Γ: rapid temperature response (dust-like)
- Typical values: 50-2500 J m⁻² K⁻¹ s⁻¹/² for a planetary surface

# Note
The unit is sometimes called "tiu" (thermal inertia unit).
"""
thermal_inertia(k, ρ, Cₚ) = @. √(k * ρ * Cₚ)


"""
    thermal_diffusivity(k, ρ, Cp) -> α

Calculate the thermal diffusivity of a material.

# Arguments
- `k::Real`  : Thermal conductivity [W/m/K]
- `ρ::Real`  : Material density [kg/m³]
- `Cₚ::Real` : Heat capacity [J/kg/K]

# Returns
- `α::Real` : Thermal diffusivity [m²/s]

# Mathematical Formula
```math
\\alpha = \\frac{k}{\\rho C_p}
```

# Physical Meaning
- Measures how quickly temperature propagates through material
- Appears in the heat diffusion equation: ∂T/∂t = α∇²T
- High α: rapid heat diffusion
- Low α: slow heat diffusion
"""
thermal_diffusivity(k, ρ, Cₚ) = @. k / (ρ * Cₚ)


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
