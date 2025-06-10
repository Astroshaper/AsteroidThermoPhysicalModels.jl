#=
AsteroidThermoPhysicalModels.jl

Main module file for the AsteroidThermoPhysicalModels package.
This package provides comprehensive thermophysical modeling capabilities for asteroids,
including temperature distribution calculations and non-gravitational effects
(Yarkovsky and YORP effects).

Key features:
- 1D heat conduction with multiple numerical solvers
- Self-shadowing and self-heating effects
- Binary asteroid support with mutual effects
- Non-gravitational force calculations
=#

module AsteroidThermoPhysicalModels

using AsteroidShapeModels
using CSV
using DataFrames
using LinearAlgebra
using ProgressMeter
using StaticArrays
using Statistics

# Physical constants used throughout the package
const SOLAR_CONST = 1366.0   # Solar constant at 1 au, Φ☉ [W/m²]
                             # Average solar irradiance at Earth's mean distance
const c₀ = 299792458.0       # Speed of light in vacuum [m/s]
                             # Exact value (defined)
const σ_SB = 5.670374419e-8  # Stefan–Boltzmann constant [W/m²/K⁴]
                             # σ = 2π⁵k⁴/(15h³c²)
const h = 6.62607015e-34     # Planck constant [J⋅s]
                             # Exact value (defined in SI since 2019)
const k_B = 1.380649e-23     # Boltzmann constant [J/K]
                             # Exact value (defined in SI since 2019)
const au2m = 149597870700    # 1 astronomical unit [m]
                             # IAU 2012 definition (exact)
const m2au = 1/au2m          # Conversion factor: meters to au

# Re-export types and functions from AsteroidShapeModels.jl
export ShapeModel, VisibleFacet
export load_shape_obj, load_shape_grid
export face_center, face_normal, face_area
export polyhedron_volume, equivalent_radius, maximum_radius, minimum_radius
export view_factor, raycast, find_visiblefacets!, isilluminated

include("thermo_params.jl")
include("TPM.jl")

# Alias for the types defined in `TPM.jl`
# `TPM` is a abbreviation for a "thermophysical model".
const AbstractAsteroidTPM     = AbstractAsteroidThermoPhysicalModel
const SingleAsteroidTPM       = SingleAsteroidThermoPhysicalModel
const BinaryAsteroidTPM       = BinaryAsteroidThermoPhysicalModel
const SingleAsteroidTPMResult = SingleAsteroidThermoPhysicalModelResult
const BinaryAsteroidTPMResult = BinaryAsteroidThermoPhysicalModelResult
export AbstractAsteroidThermoPhysicalModel, SingleAsteroidThermoPhysicalModel, BinaryAsteroidThermoPhysicalModel
export AbstractAsteroidTPM, SingleAsteroidTPM, BinaryAsteroidTPM

include("heat_conduction.jl")
include("heat_conduction_analytical.jl")
include("energy_flux.jl")
include("non_grav.jl")
include("thermal_radiation.jl")
export thermal_skin_depth, thermal_inertia, init_temperature!, run_TPM!

include("roughness.jl")

end # module AsteroidThermoPhysicalModels
