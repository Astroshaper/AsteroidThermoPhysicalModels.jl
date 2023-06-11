# See https://github.com/Astroshaper/Astroshaper-examples/tree/main/TPM_Ryugu for more information.
@testset "find_visiblefacets" begin
    ##= Shape model of Ryugu =##
    filename = "SHAPE_SFM_49k_v20180804.obj"
    url_shape = "https://data.darts.isas.jaxa.jp/pub/hayabusa2/paper/Watanabe_2019/$(filename)"
    filepath = joinpath("shape", filename)
    isfile(filepath) || Downloads.download(url_shape, filepath)

    shape = AsteroidThermoPhysicalModels.load_shape_obj(filepath; scale=1000, find_visible_facets=true)

    ##= Shape model of an icosahedron =##
    filepath = "icosahedron.obj"
    shape = AsteroidThermoPhysicalModels.load_shape_obj(filepath; scale=1, find_visible_facets=true)

    println("==== $(filepath) ====")
    println(shape)

    for (id, facet) in enumerate(shape.facets)
        println("==== Facet #$id ====")
        println(facet)
    end
end
