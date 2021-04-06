using Astroshaper
using Test

@testset "Astroshaper.jl" begin
    shapepath = "Julia_Vernazza_2018.obj"
    @test shape = setShapeModel(shapepath; scale=1000)
    @test findVisibleFaces!(shape)
end

