
################################################################
#                 Data format conversion
################################################################

"""
Vector of vector to 2-D matrix
"""
function VectorVector2Matrix(v)
    m = Matrix{eltype(v[end])}(undef, length(v), 3)
    for i in eachindex(v)
        m[i, :] .= v[i]
    end
    m
end

"""
    face2node(nodes, faces, data) -> node_data

Convert face-based data to node-based data

# Arguments
- `nodes` : 2-D matrix of nodes
- `faces` : 2-D matrix of faces
- `data`  : Face-based data

# Return
- `node_data` : Node-based data
"""
function face2node(nodes, faces, data)
    node_data = Vector{eltype(data)}(undef, size(nodes, 1))
    
    for node_index in eachindex(node_data)
        face_indices = findall(faces.==node_index)
        if length(face_indices) == 0
            node_data[node_index] = NaN
            continue
        end
        node_data[node_index] = mean(data[face_index[1]] for face_index in face_indices)
    end
    node_data
end


################################################################
#              3D visualization of polyhedral shape
################################################################

"""
    showshape(shape; data=nothing)

"""
function draw(shape::Shape; data=nothing, r̂☉=[1,0,0.])
    nodes = VectorVector2Matrix(shape.nodes)
    faces = VectorVector2Matrix(shape.faces)

    set_theme!(backgroundcolor=:black)
    
    if data == nothing
        color = :gray
    elseif data == :radius
        color = [norm(v) for v in eachrow(nodes)]
    elseif data == :temperature
        surf_temps = surface_temperature(shape)
        color = face2node(nodes, faces, surf_temps)
    elseif data == :illumination
        illumination = Float64.(isIlluminated(r̂☉, shape))
        color = face2node(nodes, faces, illumination)
    else
        color = face2node(nodes, faces, data)
    end

    scene = mesh(nodes, faces, color=color)
    display(scene)
end


function draw(shape1::Shape, shape2::Shape)
    faces1 = VectorVector2Matrix(shape1.faces)
    faces2 = VectorVector2Matrix(shape2.faces)

    nodes1 = VectorVector2Matrix(shape1.nodes)

    r = [2600, 0, 0.]  # 平行移動 + オイラー角変換
    # ϕ1 =
    # ϕ2 =
    # ϕ3 =

    nodes2 = [node + r for node in shape2.nodes]
    nodes2 = VectorVector2Matrix(nodes2)
    

    set_theme!(backgroundcolor=:black)
    println("更新反映される？")
    

    scene = poly(nodes1, faces1, color=:gray, strokecolor=:black, strokewidth=1)
    poly!(nodes2, faces2, color=:gray, strokecolor=:black, strokewidth=1)
    display(scene)
end
