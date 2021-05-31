module Astroshaper

using LinearAlgebra
using StaticArrays
using StructArrays

using GLMakie  # 3D visulaization

include("constants.jl")
export AU, G, GM☉, M☉, SOLAR_CONST, c₀, σ_SB

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

include("shape.jl")
export Shape, setShapeModel, findVisibleFaces!, showshape

include("YORP.jl")
export getNetTorque, getNetTorque_shadowing, torque2rate, getTimeScale

include("thermophysics.jl")

################################################################

include("nbody.jl")
export Particle, setParticles, setOrigin!, getBaryCenter, setOrigin2BaryCenter!

include("Hermite4.jl")
export run_Hermite4, Initialize!
export run_sim, forward!, predict!, evaluate!, evaluate_by_predictor!, evaluate_by_corrector!
export collect!, prepare!, get_Δt_Aarseth, get_Δt_initial

include("Euler.jl")
export run_Euler, getParticlesCOM, setOrigin2COM!

end # module Astroshaper
