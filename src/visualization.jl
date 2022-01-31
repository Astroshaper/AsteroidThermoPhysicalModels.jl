
################################################################
#                         Plot orbits
################################################################


function spkpos_df(targ, ets::AbstractVector, ref, abcorr, obs)
    df = DataFrame(et=Float64[], x=Float64[], y=Float64[], z=Float64[], lt=Float64[])
 
    for et in ets
        pos, lt = SPICE.spkpos(targ, et, ref, abcorr, obs)
        push!(df, (et, pos..., lt))
    end

    df
end


function plot_orbits(kernels, bodies, orbit)

    for kernel in kernels
        SPICE.furnsh(kernel)
    end

    ref    = "ECLIPJ2000"
    abcorr = "none"
    obs    = "sun"

    et_start = 0.
    et_end   = SPICE.jyear() * 2
    Δet      = SPICE.spd() * 10  # seconds per day
    ets = collect(et_start:Δet:et_end);

    df = spkpos_df("Earth", ets, ref, abcorr, obs)
    df[:, [:x, :y, :z]] .= SPICE.convrt.(df[:, [:x, :y, :z]], "km", "au")
    scene = meshscatter(df.x, df.y, df.z, markersize=0.01, color="blue")

    # for body in bodies
    #     # GM = SPICE.bodvrd(body, :GM)[1]
    #     df = spkpos_df(body, ets, "J2000", "none", "sun")
    #     df[:, [:x, :y, :z]] .= SPICE.convrt.(df[:, [:x, :y, :z]], "km", "au")
    #     meshscatter!(df.x, df.y, df.z, markersize=0.01)
    # end

    ts = collect(0:orbit.T/400:orbit.T)
    df = DataFrame(et=Float64[], x=Float64[], y=Float64[], z=Float64[])
    for t in ts
        u = solveKeplerEquation2(orbit, t)
        r = get_r(orbit, u)
        r = orbit_to_inertia(r, orbit)
        push!(df, (t, r...))
    end
    df[:, [:x, :y, :z]] .= SPICE.convrt.(df[:, [:x, :y, :z]], "m", "au")
    meshscatter!(df.x, df.y, df.z, markersize=0.01, color="orange")
    
    set_theme!(backgroundcolor=:black)
    display(scene)
end




################################################################
#                    Data format conversion
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
    

    scene = poly(nodes1, faces1, color=:gray, strokecolor=:black, strokewidth=1)
    poly!(nodes2, faces2, color=:gray, strokecolor=:black, strokewidth=1)
    display(scene)
end
