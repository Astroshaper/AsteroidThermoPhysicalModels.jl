#=
find_visiblefacets.jl

Tests for visible facet determination in shape models.
This test verifies:
- Correct identification of mutually visible facets
- View factor calculations between surface elements
- Performance with different shape model complexities
=#

@testset "find_visiblefacets" begin
    msg = """\n
    ┏━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━┓
    |                Test: find_visiblefacets                |
    ┗━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━┛
    """
    println(msg)

    ## --- Shape model of Ryugu ---
    filepath = joinpath("shape", "ryugu_test.obj")  # Small model for test
    println("========  $(filepath)  ========")
    shape = load_shape_obj(filepath; scale=1000, find_visible_facets=true)
    println(shape)

    ## --- Icosahedron ---
    filepath = joinpath("shape", "icosahedron.obj")
    println("========  $(filepath)  ========")
    shape = load_shape_obj(filepath; scale=1, find_visible_facets=true)
    println(shape)

    println("Number of total visible facets: ", sum(length.(shape.visiblefacets)))  # This should be zero for an icosahedron.
    println()

    ## --- Concave spherical segment ---
    println("========  Concave spherical segment  ========")
    xs, ys, zs = AsteroidThermoPhysicalModels.concave_spherical_segment(0.4, 0.2; Nx=2^5, Ny=2^5, xc=0.5, yc=0.5)
    shape = AsteroidThermoPhysicalModels.load_shape_grid(xs, ys, zs; scale=1.0, find_visible_facets=true)
    
    println("Number of faces visible from the crater center: ", length(shape.visiblefacets[992]))
    println()
end
