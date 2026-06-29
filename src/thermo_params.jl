#=
thermo_params.jl

Thermophysical parameter definitions and calculations.
This file contains functions for computing key thermal properties:
- Thermal skin depth
- Thermal inertia
- Thermal diffusivity
Also defines the `ThermoParams` structure for storing material properties
and discretization parameters.
=#

# ╔═══════════════════════════════════════════════════════════════════╗
# ║              Thermal skin depth & Thermal inertia                 ║
# ╚═══════════════════════════════════════════════════════════════════╝

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


# ╔═══════════════════════════════════════════════════════════════════╗
# ║               Struct for thermophysical parameters                ║
# ╚═══════════════════════════════════════════════════════════════════╝

abstract type AbstractThermoParams end

"""
    struct ThermoParams <: AbstractThermoParams

Material thermal properties per facet.

For a uniform surface, pass a single `ThermoParams` constructed from scalar values.
For a non-uniform surface, pass a `ThermoParams` constructed from `Vector{Float64}` values
of length `n_face`.

# Fields
- `conductivity`    : Thermal conductivity for each facet [W/m/K]
- `density`         : Density for each facet [kg/m³]
- `heat_capacity`   : Heat capacity for each facet [J/kg/K]
- `reflectance_vis` : Reflectance in visible light for each facet [-]
- `reflectance_ir`  : Reflectance in thermal infrared for each facet [-]
- `emissivity`      : Emissivity for each facet [-]
"""
struct ThermoParams <: AbstractThermoParams
    conductivity    ::Vector{Float64}
    density         ::Vector{Float64}
    heat_capacity   ::Vector{Float64}
    reflectance_vis ::Vector{Float64}
    reflectance_ir  ::Vector{Float64}
    emissivity      ::Vector{Float64}
end


"""
    ThermoParams(conductivity, density, heat_capacity, reflectance_vis, reflectance_ir, emissivity)

Construct `ThermoParams` from scalar values (uniform surface).

# Arguments
- `conductivity`    : Thermal conductivity [W/m/K]
- `density`         : Density [kg/m³]
- `heat_capacity`   : Heat capacity [J/kg/K]
- `reflectance_vis` : Reflectance in visible light [-]
- `reflectance_ir`  : Reflectance in thermal infrared [-]
- `emissivity`      : Emissivity [-]
"""
function ThermoParams(
    conductivity    ::Float64,
    density         ::Float64,
    heat_capacity   ::Float64,
    reflectance_vis ::Float64,
    reflectance_ir  ::Float64,
    emissivity      ::Float64,
)
    ThermoParams([conductivity], [density], [heat_capacity], [reflectance_vis], [reflectance_ir], [emissivity])
end


"""
    ThermoParams(; conductivity, density, heat_capacity, reflectance_vis, reflectance_ir, emissivity)

Construct `ThermoParams` from keyword arguments.
Accepts the same types as the positional constructors (scalar `Float64` or `Vector{Float64}`).

# Keyword Arguments
- `conductivity`    : Thermal conductivity [W/m/K]
- `density`         : Density [kg/m³]
- `heat_capacity`   : Heat capacity [J/kg/K]
- `reflectance_vis` : Reflectance in visible light [-]
- `reflectance_ir`  : Reflectance in thermal infrared [-]
- `emissivity`      : Emissivity [-]
"""
function ThermoParams(;
    conductivity,
    density,
    heat_capacity,
    reflectance_vis,
    reflectance_ir,
    emissivity,
)
    ThermoParams(conductivity, density, heat_capacity, reflectance_vis, reflectance_ir, emissivity)
end


"""
    subsolar_temperature(r☉, R_vis, ε) -> Tₛₛ

Calculate the subsolar equilibrium temperature at a given heliocentric distance.

# Arguments
- `r☉`    : Sun's position vector in the asteroid-fixed frame [m]
- `R_vis` : Visible-light reflectance (Bond albedo) [-]
- `ε`     : Emissivity [-]

# Returns
- `Tₛₛ::Float64` : Subsolar point temperature [K]

# Notes
Assumes instantaneous radiative equilibrium (zero thermal inertia). Useful as an
upper bound for surface temperatures and as an initial guess for `T₀`.

# Mathematical Formula
```math
T_{ss} = \\left[\\frac{(1 - A) \\Phi_\\odot}{\\varepsilon \\sigma}\\right]^{1/4}
```
where ``\\Phi_\\odot = \\Phi_0 / r^2`` is the solar flux at heliocentric distance ``r``.
"""
function subsolar_temperature(r☉, R_vis, ε)
    Φ = SOLAR_CONST / (norm(r☉) * m2au)^2
    return ((1 - R_vis) * Φ / (ε * σ_SB))^(1/4)
end
