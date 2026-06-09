#=
tpm_run.jl

Temperature initialization and thermophysical model simulation functions.
=#

# ╔═══════════════════════════════════════════════════════════════════╗
# ║                   Initialize temperatures                         ║
# ╚═══════════════════════════════════════════════════════════════════╝


"""
    subsolar_temperature(r☉, R_vis, ε) -> Tₛₛ
    subsolar_temperature(r☉, params::AbstractThermoParams) -> Tₛₛ

Calculate the subsolar temperature on an asteroid at a given heliocentric distance.

# Arguments
- `r☉::Vector`                   : Sun's position vector in the asteroid's fixed frame [m]
- `R_vis::Real`                  : Visible light reflectance (albedo) [-]
- `ε::Real`                      : Emissivity [-]
- `params::AbstractThermoParams` : Thermal parameters (alternative input)

# Returns
- `Tₛₛ::Float64` : Subsolar point temperature [K]

# Mathematical Formula
Assuming radiative equilibrium with zero thermal conductivity (zero thermal inertia):
```math
T_{ss} = \\left[\\frac{(1 - R_{vis}) \\Phi_\\odot}{\\varepsilon \\sigma}\\right]^{1/4}
```
where:
- ``\\Phi_\\odot = \\Phi_0 / r^2`` is the solar flux at distance r [au]
- ``\\Phi_0 = 1366`` W/m² is the solar constant at 1 au
- ``\\sigma`` is the Stefan-Boltzmann constant

# Notes
- This gives the maximum temperature for a non-rotating asteroid
- Actual subsolar temperature may be lower due to thermal inertia
- Valid for airless bodies with negligible heat conduction
"""
subsolar_temperature(r☉, params::AbstractThermoParams) = subsolar_temperature(r☉, params.reflectance_vis, params.emissivity)

function subsolar_temperature(r☉, R_vis, ε)
    Φ = SOLAR_CONST / (norm(r☉) * m2au)^2  # Energy flux at the solar distance [W/m²]
    Tₛₛ = ((1 - R_vis) * Φ / (ε * σ_SB))^(1/4)

    return Tₛₛ
end


"""
    init_temperature!(stpm::SingleAsteroidThermoPhysicalModel, T₀::Real)

Initialize all temperature cells at the given temperature `T₀`

# Arguments
- `stpm` : Thermophysical model for a single asteroid
- `T₀`   : Initial temperature of all cells [K]
"""
function init_temperature!(stpm::SingleAsteroidThermoPhysicalModel, T₀::Real)
    stpm.temperature .= T₀
end


"""
    init_temperature!(btpm::BinaryAsteroidThermoPhysicalModel, T₀::Real)

Initialize all temperature cells at the given temperature `T₀`

# Arguments
- `btpm` : Thermophysical model for a binary asteroid
- `T₀`   : Initial temperature of all cells [K]
"""
function init_temperature!(btpm::BinaryAsteroidThermoPhysicalModel, T₀::Real)
    init_temperature!(btpm.pri, T₀)
    init_temperature!(btpm.sec, T₀)
end

