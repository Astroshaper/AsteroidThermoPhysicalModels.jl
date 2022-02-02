module Astroshaper

using LinearAlgebra
using StaticArrays
using StructArrays
using Statistics

import SPICE

# using Distributed
# using SharedArrays

using DataFrames
using Parameters
using GLMakie  # 3D visulaization

using FileIO
using JLD2


include("constants.jl")
export AU, G, GM☉, M☉, SOLAR_CONST, c₀, σ_SB
export MERCURY, VENUS, EARTH, MARS, JUPITER, SATURN, URANUS, NEPTUNE
export CERES, PLUTO, ERIS
export MOON
export RYUGU

include("obj.jl")
export loadobj

include("coordinates.jl")
export rotateX, rotateY, rotateZ
export rotateX!, rotateY!, rotateZ!

include("facet.jl")
include("shape.jl")
include("visualization.jl")
export Shape
export draw

include("kepler.jl")
include("spin.jl")
export OrbitalElements
export ref_to_orb!, ref_to_orb
export orb_to_ref!, orb_to_ref
export solveKeplerEquation1, solveKeplerEquation2, u2ν, heliocentric_distance
export SpinParams

include("TPM.jl")
include("thermophysics.jl")
include("YORP.jl")
export ParamsThermo
export torque2rate, YORP_timescale, run_YORP, run_Yarkovsky

include("polyhedron.jl")

include("binary.jl")
export run_binary

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
