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

include("coordinates.jl")
export rotateX, rotateY, rotateZ
export rotateX!, rotateY!, rotateZ!

include("facet.jl")
include("shape.jl")
export ShapeModel

include("kepler.jl")
include("spin.jl")
export OrbitalElements
export ref_to_orb!, ref_to_orb
export orb_to_ref!, orb_to_ref
export solveKeplerEquation1, solveKeplerEquation2, u2ν, heliocentric_distance
export SpinParams

include("thermophysics.jl")
include("TPM.jl")
include("YORP.jl")
include("Yarkovsky.jl")
export ThermoParams
export init_temps_zero!, run_TPM!
export analyze_YORP, YORP_timescale
export run_YORP, run_Yarkovsky

include("roughness.jl")
export crater_curvature_radius, concave_spherical_segment, parallel_sinusoidal_trenches

include("polyhedron.jl")
export loadmesh, nvertices, nfaces, nedges
export getSurfaceGravity

include("binary.jl")
export MutualOrbit, Binary, run_binary_TPM!


################################################################


include("nbody.jl")
export AbstractParticle, SimpleParticle, HermiteParticle
export run_nbody!
export load_snapshot
export setParticles, addParticle!, getBaryCenter, setOrigin2BaryCenter!
export getTotalEnergy, getKineticEnergy, sumKineticEnergy, getPotentialEnergy, sumPotentialEnergy

include("Hermite4.jl")
export run_Hermite4, run_Hermite4_test, initialize!

include("Euler.jl")
export run_Euler

include("leapfrog.jl")
export run_leapfrog

end # module AsteroidThermoPhysicalModels
