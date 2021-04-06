module Astroshaper


using LinearAlgebra
using StaticArrays

include("constants.jl")
export AU, G, GM☉, SOLAR_CONST, c₀

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
export Shape, setShapeModel, findVisibleFaces!

include("YORP.jl")
export getNetTorque, getNetTorque_shadowing, torque2rate, getTimeScale


end # module Astroshaper
