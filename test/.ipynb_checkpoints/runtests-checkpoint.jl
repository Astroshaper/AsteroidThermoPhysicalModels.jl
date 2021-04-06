using Astroshaper
using Test

@testset "Astroshaper.jl" begin
    shapepath = "itokawa_v2718_f5430.obj"
    @test shape = setShapeModel(shapepath; scale=1000)
    @test findVisibleFaces!(shape)
end

