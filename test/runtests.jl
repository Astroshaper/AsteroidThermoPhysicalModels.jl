using ThermoPhysicalModeling
using Test
using Aqua
using JLD2
using SPICE
using Git
using Downloads
using Statistics
using LinearAlgebra

Aqua.test_all(ThermoPhysicalModeling, ambiguities=false)

include("TPM_Ryugu.jl")
include("TPM_Didymos.jl")
include("Keplerian_motion.jl")
