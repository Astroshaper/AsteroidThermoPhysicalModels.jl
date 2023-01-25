using Astroshaper
using Test
using Aqua
using JLD2
using SPICE
using Downloads
using Statistics

Aqua.test_all(Astroshaper, ambiguities=false)

@testset "Astroshaper.jl" begin
    # shapepath = "./ryugu_test.obj"
    # shape = Shape(shapepath; scale=1000, find_visible_facets=true, save_shape=false)
    # @test typeof(shape) == Shape
    @test true
end

include("TPM_Ryugu.jl")
