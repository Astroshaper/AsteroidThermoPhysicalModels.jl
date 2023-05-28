

function isobj(filepath)
    base, ext = splitext(filepath)
    return ext == ".obj"
end


"""
    loadobj(shapepath::String; scale=1, message=true) -> nodes, faces
"""
function loadobj(shapepath::String; scale=1, static=true, message=true)

    nodes = SVector{3, Float64}[]
    faces = SVector{3, Int64}[]

    open(shapepath, "r") do f
        for line in eachline(f)
            line == "" && continue  # Skip blank lines
            data = split(line)

            if data[1] == "v"
                node = parse.(Float64, data[2:4])  # (x, y, z) coordinates [km]
                push!(nodes, node)
            elseif data[1] == "f"
                face = parse.(Int64, data[2:4])  # indices of three vertices
                push!(faces, face)
            end
        end
    end

    nodes *= scale  # if scale is 1000, converted [km] to [m]
    
    if static == true
        nodes = [SVector{3,Float64}(node) for node in nodes]
        faces = [SVector{3,Int64}(face) for face in faces]
    end
    
    if message == true
        println("+-----------------------------+")
        println("|        Load OBJ file        |")
        println("+-----------------------------+")
        println(" Nodes: ", length(nodes))
        println(" Faces: ", length(faces))
    end

    return nodes, faces
end
