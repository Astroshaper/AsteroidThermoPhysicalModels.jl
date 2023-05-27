

"""
    ShapeModel{T}

A polyhedral shape model of an asteroid.

# Fields
- `num_node`   : Number of nodes
- `num_face`   : Number of faces
- `nodes`      : 1-D array of node positions
- `faces`      : 1-D array of vertex indices of faces

- `facets`     : 1-D array of surface facets (`Facet`)

- `VOLUME`     : Volume
- `RADIUS_EQ`  : Equivalent radius of a sphere with the same volume
- `RADIUS_MAX` : Maximum radius
- `RADIUS_MIN` : Minimum radius
- `COF`        : Center-of-figure

- `force`      : Thermal recoil force at body-fixed frame (Yarkovsky effect)
- `torque`     : Thermal recoil torque at body-fixed frame (YORP effect)
"""
struct ShapeModel{T}
    num_node  ::Int
    num_face  ::Int
    nodes     ::Vector{SVector{3, Float64}}
    faces     ::Vector{SVector{3, Int}}
    facets    ::T
    VOLUME    ::Float64
    RADIUS_EQ ::Float64
    RADIUS_MAX::Float64
    RADIUS_MIN::Float64
    force     ::MVector{3, Float64}
    torque    ::MVector{3, Float64}
end


function Base.show(io::IO, shape::ShapeModel)
    msg = "Shape model\n"
    msg *= "-----------\n"
    msg *= "Nodes             : $(shape.num_node)\n"
    msg *= "Faces             : $(shape.num_face)\n"
    msg *= "Volume            : $(shape.VOLUME)\n"
    msg *= "Equivalent radius : $(shape.RADIUS_EQ)\n"
    msg *= "Maximum radius    : $(shape.RADIUS_MAX)\n"
    msg *= "Minimum radius    : $(shape.RADIUS_MIN)\n"
    print(io, msg)
end


function ShapeModel(shapepath; scale=1, find_visible_facets=false, save_shape=false)
    # TODO: use MeshIO.jl
    
    ext = splitext(shapepath)[2]
    
    if ext == ".obj"
        nodes, faces = loadobj(shapepath; scale=scale, static=true, message=false)
        
        num_node = length(nodes)
        num_face = length(faces)
        
        facets = getfacets(nodes, faces)
        find_visible_facets && findVisibleFacets!(facets)
        
        VOLUME     = getvolume(facets)
        RADIUS_EQ  = equivalent_radius(VOLUME)
        RADIUS_MAX = maximum_radius(nodes)
        RADIUS_MIN = minimum_radius(nodes)
        
        force  = zero(MVector{3, Float64})
        torque = zero(MVector{3, Float64})
        
        shape = ShapeModel(
            num_node, num_face, nodes, faces, facets,
            VOLUME, RADIUS_EQ, RADIUS_MAX, RADIUS_MIN, force, torque
        )
        save_shape && save(splitext(shapepath)[1] * ".jld2", Dict("shape" => shape))

    elseif ext == ".jld2"
        shape = load(shapepath)["shape"]
    else
        return ArgumentError("Give a filepath of *.obj or *.jld2.")
    end

    return shape
end

equivalent_radius(VOLUME) = (3VOLUME/4π)^(1/3)
equivalent_radius(shape::ShapeModel) = equivalent_radius(shape.VOLUME)

maximum_radius(nodes) = maximum(norm.(nodes))
maximum_radius(shape::ShapeModel) = maximum_radius(shape.nodes)

minimum_radius(nodes) = minimum(norm.(nodes))
minimum_radius(shape::ShapeModel) = minimum_radius(shape.nodes)

findVisibleFacets!(shape::ShapeModel) = findVisibleFacets!(shape.facets)
isIlluminated(obs::Facet, r̂☉, shape::ShapeModel) = isIlluminated(obs, r̂☉, shape.facets)
isIlluminated(r̂☉, shape::ShapeModel) = [isIlluminated(obs, r̂☉, shape) for obs in shape.facets]

surface_temperature(shape::ShapeModel) = [facet.temps[begin] for facet in shape.facets]


################################################################
#                      Shape properites
################################################################


"""
    getvolume(facets) -> VOLUME

Calculate volume of a polyhedral
"""
getvolume(facets) = sum(((facet.A × facet.B) ⋅ facet.C) / 6 for facet in facets)
