using Astroshaper
using Test

@testset "Astroshaper.jl" begin
    shapepath = "./ryugu_test.obj"
    shape = setShapeModel(shapepath; scale=1000)
    @time Astroshaper.sumTorqueOverSurface(shape, 1200, [1,0,0.])
    @test true
end

