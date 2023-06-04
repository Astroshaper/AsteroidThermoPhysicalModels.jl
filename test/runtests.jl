using AsteroidThermoPhysicalModels
using Test
using Aqua
using JLD2
using SPICE
using Git
using Downloads
using Statistics
using LinearAlgebra
using StaticArrays
using Rotations

ENABLE_JLD = true
Aqua.test_all(AsteroidThermoPhysicalModels, ambiguities=false)

include("TPM_Ryugu.jl")
include("TPM_Didymos.jl")
include("thermal_property_distribution.jl")
