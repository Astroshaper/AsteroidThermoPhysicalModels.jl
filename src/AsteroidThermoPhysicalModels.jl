module AsteroidThermoPhysicalModels

using LinearAlgebra
using StaticArrays
using Statistics

import SPICE

using DataFrames
using ProgressMeter

using FileIO
using CSV

import MeshIO
import GeometryBasics

const SOLAR_CONST = 1366.0   # Solar constant, Φ☉ [W/m^2]
const c₀ = 299792458.0       # Speed of light [m/s]
const σ_SB = 5.670374419e-8  # Stefan–Boltzmann constant [W/m^2/K^4]
const h = 6.62607015e-34     # Planck constant [J⋅s]
const k_B = 1.380649e-23     # Boltzmann constant [J/K]

include("obj.jl")
include("shape.jl")
include("facet.jl")
export ShapeModel, load_shape_obj

include("thermo_params.jl")
include("TPM.jl")
include("heat_conduction.jl")
include("energy_flux.jl")
include("non_grav.jl")
include("radiation.jl")
export thermal_skin_depth, thermal_inertia, init_temperature!, run_TPM!

include("roughness.jl")

end # module AsteroidThermoPhysicalModels
