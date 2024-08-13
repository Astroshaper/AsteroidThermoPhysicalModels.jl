

function isobj(filepath)
    base, ext = splitext(filepath)
    return ext == ".obj"
end


"""
    loadobj(shapepath::String; scale=1, message=true) -> nodes, faces
"""
function loadobj(shapepath::String; scale = 1, message = true)

    nodes = SVector{3,Float64}[]
    faces = SVector{3,Int64}[]

    mesh = load(shapepath)
    nodes = Vector{SVector{3, Float64}}(GeometryBasics.coordinates(mesh))
    faces = [SVector{3,Int}(convert.(Int, face)) for face in GeometryBasics.faces(mesh)]

    nodes *= scale  # if scale is 1000, converted [km] to [m]

    if message == true
        println("+-----------------------------+")
        println("|        Load OBJ file        |")
        println("+-----------------------------+")
        println(" Nodes: ", length(nodes))
        println(" Faces: ", length(faces))
    end

    return nodes, faces
end
