module Astroshaper

using LinearAlgebra
using StaticArrays
using StructArrays

using DataFrames
using Parameters
using GLMakie  # 3D visulaization

using FileIO
using JLD2


include("constants.jl")
export AU, G, GM☉, M☉, SOLAR_CONST, c₀, σ_SB
export M_Mercury, GM_Mercury
export M_Venus,   GM_Venus
export M_Earth,   GM_Earth
export M_Moon,    GM_Moon
export M_Mars,    GM_Mars
export M_Ceres,   GM_Ceres
export M_Jupiter, GM_Jupiter
export M_Saturn,  GM_Saturn
export M_Uranus,  GM_Uranus
export M_Neptune, GM_Neptune
export M_Pluto,   GM_Pluto
export M_Eris,    GM_Eris

include("obj.jl")
export loadobj

include("coordinates.jl")
export rotateX, rotateY, rotateZ
export rotateX!, rotateY!, rotateZ!

include("kepler.jl")
export OrbitalElements
export ref_to_orb!, ref_to_orb
export orb_to_ref!, orb_to_ref

include("spin.jl")
export Spin, setSpinParams

include("smesh.jl")
export SMesh

include("Facet.jl") ## <-------------------------------------

include("shape.jl")
export Shape, setShapeModel, findVisibleFaces!, showshape, equivalent_radius, get_surface_temperature

include("YORP.jl")
export getNetTorque, getNetTorque_shadowing, torque2rate, getTimeScale, run_YORP, run_Yarkovsky

include("thermophysics.jl")
export ParamsThermo

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
