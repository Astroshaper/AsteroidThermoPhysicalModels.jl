module AsteroidThermoPhysicalModels

using LinearAlgebra
using StaticArrays
using StructArrays
using Statistics

import SPICE

# using Distributed
# using SharedArrays

using DataFrames
using Parameters
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
export loadobj

include("facet.jl")
include("shape.jl")
export ShapeModel

include("thermophysics.jl")
include("TPM.jl")
include("non_grav.jl")
export ThermoParams
export init_temps_zero!, run_TPM!

include("roughness.jl")
export crater_curvature_radius, concave_spherical_segment, parallel_sinusoidal_trenches

include("polyhedron.jl")
export loadmesh, nvertices, nfaces, nedges
export getSurfaceGravity

end # module AsteroidThermoPhysicalModels
