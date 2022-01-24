module Astroshaper

using LinearAlgebra
using StaticArrays
using StructArrays
using Statistics

# using Distributed
# using SharedArrays

using DataFrames
using Parameters
using GLMakie  # 3D visulaization

using FileIO
using JLD2


include("constants.jl")
export AU, G, GM☉, M☉, SOLAR_CONST, c₀, σ_SB
export MERCURY, VENUS, EARTH, MOON, MARS, CERES, JUPITER, SATURN, URANUS, NEPTUNE, PLUTO, ERIS
export RYUGU

include("obj.jl")
export loadobj

include("coordinates.jl")
export rotateX, rotateY, rotateZ
export rotateX!, rotateY!, rotateZ!

include("kepler.jl")
export OrbitalElements
export ref_to_orb!, ref_to_orb
export orb_to_ref!, orb_to_ref
export solveKeplerEquation1, solveKeplerEquation2, u2ν, heliocentric_distance

include("spin.jl")
export Spin, setSpinParams

include("facet.jl")

include("shape.jl")
export Shape, draw

include("YORP.jl")
export getNetTorque, getNetTorque_shadowing, torque2rate, getTimeScale, run_YORP, run_Yarkovsky

include("thermophysics.jl")
export ParamsThermo

include("polyhedron.jl")

################################################################

include("nbody.jl")
export AbstractParticle, SimpleParticle, HermiteParticle
export run_nbody!
export load_snapshot
export setParticles, addParticle!, setOrigin!, getBaryCenter, setOrigin2BaryCenter!
export getTotalEnergy, getKineticEnergy, sumKineticEnergy, getPotentialEnergy, sumPotentialEnergy

include("Hermite4.jl")
export run_Hermite4, run_Hermite4_test, initialize!

include("Euler.jl")
export run_Euler

include("leapfrog.jl")
export run_leapfrog

# include("RungeKutta.jl")

end # module Astroshaper
