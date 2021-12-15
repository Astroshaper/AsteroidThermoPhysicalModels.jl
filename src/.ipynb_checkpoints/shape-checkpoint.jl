

"""
    Shape{T1, T2, T3, T4, T5, T6, T7}

A polyhedral shape model of an asteroid.

# Fields

- `num_node` : Number of nodes
- `num_face` : Number of faces
- `nodes`    : 1-D array of node positions
- `faces`    : 1-D array of vertex indices of faces

- `facets`   : 1-D array of surface facets (`Facet`)

- `AREA`     : Surface area
- `VOLUME`   : Volume
- `COF`      : Center-of-figure
- `MOI`      : Moment of inertia tensor

– `force`    : Thermal recoil force at body-fixed frame (Yarkovsky effect)
- `torque`   : Thermal recoil torque at body-fixed frame (YORP effect)
- `Tz⁺`      : Pre-allocated vector for update of temperature profile on each facets (`Facet`)
"""
struct Shape{T1, T2, T3, T4, T5, T6, T7, T8, T9, T10}
    num_node::T1
    num_face::T1
    nodes   ::T2
    faces   ::T3

    facets::T4

    AREA  ::T5
    VOLUME::T5
    COF   ::T6
    MOI   ::T7

    force ::T8
    torque::T9
    Tz⁺   ::T10
end


function Base.show(io::IO, shape::Shape)
    println(io, "Shape model")
    println("-----------")

    println("Nodes            : ", shape.num_node)
    println("Faces            : ", shape.num_face)
    println("Surface area     : ", shape.AREA)
    println("Volume           : ", shape.VOLUME)
    println("Center-of-Figure : ", shape.COF)
    println("Inertia tensor   : ")
    println("    | Ixx Ixy Ixz |   ", shape.MOI[1, :])
    println("    | Iyx Iyy Iyz | = ", shape.MOI[2, :])
    println("    | Izx Izy Izz |   ", shape.MOI[3, :])
end


function setShapeModel(shapepath; scale=1, find_visible_facets=false, save_shape=false)
    
    ext = splitext(shapepath)[2]
    
    if ext == ".obj"
        nodes, faces = loadobj(shapepath; scale=scale, static=true, message=false)
        
        num_node = length(nodes)
        num_face = length(faces)
        
        facets = getfacets(nodes, faces)
        find_visible_facets && findVisibleFacets!(facets)
        
        AREA   = sum(facets.area)
        VOLUME = getvolume(facets)
        COF    = getCOF(facets)
        MOI    = getMOI(facets)
        
        force  = zeros(3)
        torque = zeros(3)
        Tz⁺    = similar(facets[begin].Tz)
        
        shape = Shape(
                num_node, num_face, nodes, faces, facets,
                AREA, VOLUME, COF, MOI,
                force, torque, Tz⁺
        )
        save_shape && save(splitext(shapepath)[1] * ".jld2", Dict("shape" => shape))

    elseif ext == ".jld2"
        shape = load(shapepath)["shape"]
    else
        return ArgumentError("Give a filepath of *.obj or *.jld2.")
    end

    return shape
end


getFaceCenters(shape::Shape) = getcenters(shape.smeshes)
getFaceNormals(shape::Shape) = getnormals(shape.smeshes)
getFaceAreas(shape::Shape) = getareas(shape.smeshes)

equivalent_radius(VOLUME) = (3VOLUME/4π)^(1/3)
equivalent_radius(shape::Shape) = equivalent_radius(shape.VOLUME)

findVisibleFacets!(shape::Shape) = findVisibleFacets!(shape.facets)
isIlluminated(obs::Facet, r̂☉, shape::Shape) = isIlluminated(obs, r̂☉, shape.facets)
isIlluminated(r̂☉, shape::Shape) = [isIlluminated(obs, r̂☉, shape) for obs in shape.facets]


################################################################
#                      Shape properites
################################################################


"""
    getvolume(facets) -> VOLUME

Calculate volume of a polyhedral
"""
getvolume(facets) = sum(((facet.A × facet.B) ⋅ facet.C) / 6 for facet in facets)

"""
    getCOF(facets) -> COF

Calculate center-of-figure position of a polyhedral
"""
function getCOF(facets)
    VOLUME = getvolume(facets)
    COF = zeros(3)

    for facet in facets
        @unpack A, B, C = facet
        volume = ((A × B) ⋅ C) / 6  # volume of pyramid element O-A-B-C
        center = (A + B + C) / 4    # center of pyramid element O-A-B-C
        COF += volume * center
    end
    COF / VOLUME
end

"""
    getMOI(facets) -> MOI

Calculate moment of inertia tensor of a polyhedron
"""
function getMOI(facets)
    MOI = zeros(3, 3)

    for facet in facets
        # v1, v2, v3 = m.vs
        
        # I[3, 3] += (v1[1]*v1[1] + v1[1]*v2[1] + v2[1]*v2[1] + v1[1]*v3[1] + v2[1]*v3[1] + v3[1]*v3[1] + v1[2]*v1[2] + v1[2]*v2[2] + v2[2]*v2[2] + v1[2]*v3[2] + v2[2]*v3[2] + v3[2]*v3[2]) / 60
    end
    MOI
end


################################################################
#              3D visualization of polyhedral shape
################################################################

"""
Vector of vector to 2D Matrix
"""
function VectorVector2Matrix(v)
    m = Matrix{eltype(v[end])}(undef, length(v), 3)
    for i in eachindex(v)
        m[i, :] .= v[i]
    end
    m
end


function showshape(shape)
    nodes = VectorVector2Matrix(shape.nodes)
    faces = VectorVector2Matrix(shape.faces)
    
    colors = :gray
    # colors = isIlluminated([1,0,0], shape)  # 色付けは、nodeベースしかできない？
    
    # scene = mesh(nodes, faces, color = colors, shading = false)
    set_theme!(backgroundcolor = :black)
    scene = mesh(nodes, faces, color=colors)
    display(scene)
end


get_surface_temperature(shape) = [facet.Tz[begin] for facet in shape.facets]

