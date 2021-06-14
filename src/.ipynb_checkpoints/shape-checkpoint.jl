

"""
    Shape{T1, T2, T3, T4, T5, T6, T7}

A polyhedral shape model of an asteroid.

# Fields

- `num_node` : Number of nodes
- `num_face` : Number of faces
- `nodes`    : 1-D array of node positions
- `faces`    : 1-D array of vertex indices of faces

- `smeshes`  : 1-D array of surface meshes (`SMesh`)

- `AREA`     : Surface area
- `VOLUME`   : Volume
- `COF`      : Center-of-figure
- `I`        : Moment of inertia tensor

- `τ`        : Thermal recoil torque at body-fixed frame
"""
struct Shape{T1, T2, T3, T4, T5, T6, T7, T8}
    num_node::T1
    num_face::T1
    nodes::T2
    faces::T3

    smeshes::T4

    AREA::T5
    VOLUME::T5
    COF::T6
    I::T7

    τ::T8
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
    println("    | Ixx Ixy Ixz |   ", shape.I[1, :])
    println("    | Iyx Iyy Iyz | = ", shape.I[2, :])
    println("    | Izx Izy Izz |   ", shape.I[3, :])
end


function setShapeModel(shapepath::AbstractString; scale=1, find_visible_faces=false)
    nodes, faces = loadobj(shapepath; scale=scale, static=true, message=false)

    num_node = length(nodes)
    num_face = length(faces)
    
    smeshes = getmeshes(nodes, faces)
    find_visible_faces && findVisibleFaces!(smeshes)

    AREA = sum(getareas(smeshes))
    VOLUME = getvolume(smeshes)
    COF = getCOF(smeshes)
    I = getMOI(smeshes)
    
    τ = zeros(3)

    Shape(num_node, num_face, nodes, faces, smeshes, AREA, VOLUME, COF, I, τ)
end


getFaceCenters(shape::Shape) = getcenters(shape.smeshes)
getFaceNormals(shape::Shape) = getnormals(shape.smeshes)
getFaceAreas(shape::Shape) = getareas(shape.smeshes)

findVisibleFaces!(shape::Shape) = findVisibleFaces!(shape.smeshes)
isIlluminated(obs::SMesh, r̂☉, shape::Shape) = isIlluminated(obs, r̂☉, shape.smeshes)


################################################################
#                      Shape properites
################################################################


"""
Calculate volume of a polyhedral
"""
getvolume(smeshes) = sum((((m.A × m.B) ⋅ m.C) / 6 for m in smeshes))


"""
Calculate center-of-figure of a polyhedral
"""
function getCOF(smeshes)
    VOLUME = getvolume(smeshes)
    COF = zeros(3)

    for m in smeshes
        volume = ((m.A × m.B) ⋅ m.C) / 6  # volume of pyramid element O-A-B-C
        center = (m.A + m.B + m.C) / 4    # center of pyramid element O-A-B-C
        COF += volume * center
    end
    COF / VOLUME
end


"""
Calculate moment of inertia tensor of a polyhedron
"""
function getMOI(smeshes)
    I = zeros(3, 3)

    for m in smeshes
        # v1, v2, v3 = m.vs
        
        # I[3, 3] += (v1[1]*v1[1] + v1[1]*v2[1] + v2[1]*v2[1] + v1[1]*v3[1] + v2[1]*v3[1] + v3[1]*v3[1] + v1[2]*v1[2] + v1[2]*v2[2] + v2[2]*v2[2] + v1[2]*v3[2] + v2[2]*v3[2] + v3[2]*v3[2]) / 60

    end
    I
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
    
    # scene = mesh(nodes, faces, color = colors, shading = false)
    set_theme!(backgroundcolor = :black)
    scene = mesh(nodes, faces, color=colors)
end
