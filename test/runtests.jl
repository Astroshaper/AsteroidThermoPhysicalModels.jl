#=
runtests.jl

Main test suite for `AsteroidThermoPhysicalModels.jl` package.
This file orchestrates all tests including:
- Code quality checks using Aqua.jl
- Thermophysical model tests for various asteroids
- Unit tests for individual components
- Validation against analytical solutions
=#

using AsteroidThermoPhysicalModels
using Test
using Aqua
using SPICE
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
include("test_boundary_conditions.jl")
