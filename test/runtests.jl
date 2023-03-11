using AsteroidThermoPhysicalModels
using Test
using Aqua
using JLD2
using SPICE
using Git
using Downloads
using Statistics
using LinearAlgebra

Aqua.test_all(AsteroidThermoPhysicalModels, ambiguities=false)

if !Sys.iswindows()
    # TODO: add test for Windows
    # `wget` command does not work for Windows, so the testset is disabled for now.
    include("TPM_Ryugu.jl")
end
include("TPM_Didymos.jl")
include("Keplerian_motion.jl")