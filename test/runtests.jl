using AsteroidThermoPhysicalModels
using Test
using Aqua
using SPICE
using Git
using Downloads
using Statistics
using LinearAlgebra
using StaticArrays
using Rotations
using DataFrames

Aqua.test_all(AsteroidThermoPhysicalModels, ambiguities=false)

include("find_visiblefacets.jl")
include("TPM_Ryugu.jl")
include("non-uniform_thermoparams.jl")
include("TPM_Didymos.jl")
include("heat_conduction_1D.jl")
