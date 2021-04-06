using Astroshaper
using Test

@testset "Astroshaper.jl" begin
    shapepath = "./itokawa_v2718_f5430.obj"
    shape = setShapeModel(shapepath; scale=1000)
    findVisibleFaces!(shape)
    
    @test true
end

