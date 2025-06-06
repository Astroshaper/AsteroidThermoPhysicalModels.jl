module AsteroidThermoPhysicalModels

using CSV
using DataFrames
using FileIO
using LinearAlgebra
using ProgressMeter
using StaticArrays
using Statistics

import GeometryBasics
import MeshIO

const SOLAR_CONST = 1366.0   # Solar constant, Φ☉ [W/m^2]
const c₀ = 299792458.0       # Speed of light [m/s]
const σ_SB = 5.670374419e-8  # Stefan–Boltzmann constant [W/m^2/K^4]
const h = 6.62607015e-34     # Planck constant [J⋅s]
const k_B = 1.380649e-23     # Boltzmann constant [J/K]
const au2m = 149597870700    # 1 astronomical unit is $au2m meters, same as `SPICE.convrt(1, "au", "m")`
const m2au = 1/au2m  # 1 meter is $m2au astronomical unit, same as `SPICE.convrt(1, "m", "au")`

include("obj.jl")
include("shape.jl")
include("facet.jl")
export ShapeModel, load_shape_obj

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
