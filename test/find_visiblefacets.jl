# See https://github.com/Astroshaper/Astroshaper-examples/tree/main/TPM_Ryugu for more information.
@testset "find_visiblefacets" begin
    ##= Shape model of Ryugu =##
    filename = "SHAPE_SFM_49k_v20180804.obj"
    url_shape = "https://data.darts.isas.jaxa.jp/pub/hayabusa2/paper/Watanabe_2019/$(filename)"
    filepath = joinpath("shape", filename)
    isfile(filepath) || Downloads.download(url_shape, filepath)

    shape = AsteroidThermoPhysicalModels.load_shape_obj(filepath; scale=1000, find_visible_facets=true)

    ##= Icosahedron =##
    filepath = "icosahedron.obj"
    shape = AsteroidThermoPhysicalModels.load_shape_obj(filepath; scale=1, find_visible_facets=true)

    println("==== $(filepath) ====")
    println(shape)

    total_visiblefacets = sum(length(facet.visiblefacets) for facet in shape.facets)
    println("Number of total visible facets: $total_visiblefacets")  # This should be zero for an icosahedron.

    ##= Concave spherical segment =##
    xs, ys, zs = AsteroidThermoPhysicalModels.concave_spherical_segment(0.4, 0.2; Nx=2^5, Ny=2^5, xc=0.5, yc=0.5)
    shape = AsteroidThermoPhysicalModels.load_shape_grid(xs, ys, zs; scale=1.0, find_visible_facets=true)

    println("==== Concave spherical segment ====")
    println("A facet around the crater center:")
    println(shape.facets[992])
end
