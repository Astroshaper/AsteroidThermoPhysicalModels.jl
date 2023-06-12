

"""
    ShapeModel

A polyhedral shape model of an asteroid.

# Fields
- `nodes`      : 1-D array of node positions
- `faces`      : 1-D array of vertex indices of faces
- `facets`     : 1-D array of surface facets (`Facet`)
- `force`      : Thermal recoil force at body-fixed frame (Yarkovsky effect)
- `torque`     : Thermal recoil torque at body-fixed frame (YORP effect)
"""
struct ShapeModel
    nodes     ::Vector{SVector{3, Float64}}
    faces     ::Vector{SVector{3, Int}}
    facets    ::Vector{AsteroidThermoPhysicalModels.Facet}
    force     ::MVector{3, Float64}
    torque    ::MVector{3, Float64}
end


function Base.show(io::IO, shape::ShapeModel)
    msg = "Shape model\n"
    msg *= "-----------\n"
    msg *= "Number of nodes   : $(length(shape.nodes))\n"
    msg *= "Number of faces   : $(length(shape.faces))\n"
    msg *= "Volume            : $(getvolume(shape))\n"
    msg *= "Equivalent radius : $(equivalent_radius(shape))\n"
    msg *= "Maximum radius    : $(maximum_radius(shape))\n"
    msg *= "Minimum radius    : $(minimum_radius(shape))\n"
    print(io, msg)
end

function load_shape_obj(shapepath; scale=1.0, find_visible_facets=false)
    # TODO: use MeshIO.jl
    nodes, faces = loadobj(shapepath; scale=scale, static=true, message=false)
    facets = getfacets(nodes, faces)
    find_visible_facets && find_visiblefacets!(facets)
    force  = zero(MVector{3, Float64})
    torque = zero(MVector{3, Float64})
    shape = ShapeModel(nodes, faces, facets, force, torque)
    return shape
end

function load_shape_jld(shapepath)
    shape = load(shapepath, "shape")
    return shape
end

function save_shape_jld(shapepath, shape)
    save(splitext(shapepath)[1] * ".jld2", Dict("shape" => shape))
end


"""
    load_shape_grid(xs, ys, zs; scale=1.0, find_visible_facets=false) -> shape

Convert a regular grid (x, y) to a shape model

# Arguments
- `xs::AbstractVector` : x-coordinates of grid points
- `ys::AbstractVector` : y-coordinates of grid points
- `zs::AbstractMatrix` : z-coordinates of grid points
"""
function load_shape_grid(xs::AbstractVector, ys::AbstractVector, zs::AbstractMatrix; scale=1.0, find_visible_facets=false)
    nodes, faces, facets = grid_to_facets(xs, ys, zs)
    find_visible_facets && find_visiblefacets!(facets)
    force  = zero(MVector{3, Float64})
    torque = zero(MVector{3, Float64})
    shape = ShapeModel(nodes, faces, facets, force, torque)
    return shape
end


################################################################
#                      Shape properites
################################################################

equivalent_radius(VOLUME::Real) = (3VOLUME/4π)^(1/3)
equivalent_radius(shape::ShapeModel) = equivalent_radius(getvolume(shape))

maximum_radius(nodes::Vector{<:StaticVector{3}}) = maximum(norm, nodes)
maximum_radius(shape::ShapeModel) = maximum_radius(shape.nodes)

minimum_radius(nodes) = minimum(norm.(nodes))
minimum_radius(shape::ShapeModel) = minimum_radius(shape.nodes)

find_visiblefacets!(shape::ShapeModel) = find_visiblefacets!(shape.facets)
isIlluminated(obs::Facet, r̂☉, shape::ShapeModel) = isIlluminated(obs, r̂☉, shape.facets)
isIlluminated(r̂☉, shape::ShapeModel) = [isIlluminated(obs, r̂☉, shape) for obs in shape.facets]

surface_temperature(shape::ShapeModel) = [facet.temps[begin] for facet in shape.facets]

"""
    getvolume(facets) -> VOLUME

Calculate volume of a polyhedral
"""
getvolume(facets) = sum(((facet.A × facet.B) ⋅ facet.C) / 6 for facet in facets)
getvolume(shape::ShapeModel) = getvolume(shape.facets)
