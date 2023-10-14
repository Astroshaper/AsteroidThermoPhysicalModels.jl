module AsteroidThermoPhysicalModels

using LinearAlgebra
using StaticArrays
using Statistics

import SPICE

using DataFrames
using ProgressMeter

using FileIO
using JLD2


include("constants.jl")
export AU, G, GM☉, M☉, SOLAR_CONST, c₀, σ_SB
export MERCURY, VENUS, EARTH, MARS, JUPITER, SATURN, URANUS, NEPTUNE
export CERES, PLUTO, ERIS
export MOON
export RYUGU, DIDYMOS, DIMORPHOS

include("obj.jl")
include("shape.jl")
include("facet.jl")
export ShapeModel, load_shape_obj

include("thermo_params.jl")
include("TPM.jl")
include("heat_conduction.jl")
include("energy_flux.jl")
include("non_grav.jl")
export thermal_skin_depth, thermal_inertia, init_temperature!, run_TPM!

include("roughness.jl")

end # module AsteroidThermoPhysicalModels
