# See https://github.com/Astroshaper/Astroshaper-examples/tree/main/TPM_Ryugu for more information.
@testset "find_visiblefacets" begin
    msg = """
    
    ⋅----------------------------------------⋅
    |        Test: find_visiblefacets        |
    ⋅----------------------------------------⋅
    """
    println(msg)

    ##= Shape model of Ryugu =##
    filename = "SHAPE_SFM_49k_v20180804.obj"
    url_shape = "https://data.darts.isas.jaxa.jp/pub/hayabusa2/paper/Watanabe_2019/$(filename)"
    filepath = joinpath("shape", filename)
    isfile(filepath) || Downloads.download(url_shape, filepath)

    shape = AsteroidThermoPhysicalModels.load_shape_obj(filepath; scale=1000, find_visible_facets=true)
    println("========  $(filepath)  ========")
    println(shape)

    ##= Icosahedron =##
    filepath = joinpath("shape", "icosahedron.obj")
    shape = AsteroidThermoPhysicalModels.load_shape_obj(filepath; scale=1, find_visible_facets=true)
    println("========  $(filepath)  ========")
    println(shape)

    println("Number of total visible facets: ", sum(length.(shape.visiblefacets)))  # This should be zero for an icosahedron.
    println()

    ##= Concave spherical segment =##
    println("========  Concave spherical segment  ========")
    xs, ys, zs = AsteroidThermoPhysicalModels.concave_spherical_segment(0.4, 0.2; Nx=2^5, Ny=2^5, xc=0.5, yc=0.5)
    shape = AsteroidThermoPhysicalModels.load_shape_grid(xs, ys, zs; scale=1.0, find_visible_facets=true)
    
    println(length(shape.visiblefacets[992]), " faces are visible from the crater center:")
    println()
end
