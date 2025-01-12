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

include("TPM_Ryugu/TPM_Ryugu.jl")
include("TPM_non-uniform_thermoparams/TPM_non-uniform_thermoparams.jl")
include("TPM_zero-conductivity/TPM_zero-conductivity.jl")
include("TPM_Didymos/TPM_Didymos.jl")

include("thermal_radiation.jl")
include("find_visiblefacets.jl")
include("heat_conduction_1D.jl")
