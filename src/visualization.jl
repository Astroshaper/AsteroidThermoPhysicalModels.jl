
################################################################
#                         Plot orbits
################################################################


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
    draw(shape; data=nothing)

"""
function draw(shape::ShapeModel; data=nothing, r̂☉=[1,0,0.], colormap=:viridis, strokecolor=:gray20)
    nodes = VectorVector2Matrix(shape.nodes)
    faces = VectorVector2Matrix(shape.faces)

    set_theme!(backgroundcolor=:black)
    # set_theme!(backgroundcolor=:white)
    
    if data === nothing
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

    scene = poly(nodes, faces,
        color=color, colormap=colormap,
        strokecolor=strokecolor, strokewidth=1,
        size=(1500,1500)
    )

    display(scene)

    # fig = Figure()
    # ax = Axis3(fig[1, 1], aspect=:data)

    # poly!(ax, nodes, faces,
    #     color=color, colormap=colormap,
    #     strokecolor=strokecolor, strokewidth=1,
    #     size=(1500,1500)
    # )

    # # Colorbar(fig[1, 2], limits = (minimum(surf_temps), maximum(surf_temps)), colormap=colormap)

    # display(fig)
end


"""
適度に回転させて表示する
（二重小惑星の重心を原点とする座標系）
"""
function draw(binary::Binary; color=:gray)
    @unpack shape1, shape2, mutual_orbit, spin1, spin2 = binary

    faces1 = VectorVector2Matrix(shape1.faces)
    faces2 = VectorVector2Matrix(shape2.faces)

    nodes1 = collect.(shape1.nodes)
    nodes2 = collect.(shape2.nodes)

    for (node1, node1_static) in zip(nodes1, shape1.nodes)
        node1 .= rotateZ(node1_static, -spin1.ϕ)
        node1 .+= mutual_orbit.r₁
    end

    for (node2, node2_static) in zip(nodes2, shape2.nodes)
        node2 .= rotateZ(node2_static, -spin2.ϕ)
        node2 .+= mutual_orbit.r₂
    end

    nodes1 = VectorVector2Matrix(nodes1)
    nodes2 = VectorVector2Matrix(nodes2)
    
    set_theme!(backgroundcolor=:black)

    if color == :gray
        color1 = :gray
        color2 = :gray
    elseif color == :temperature
        color1 = face2node(nodes1, faces1, surface_temperature(shape1))
        color2 = face2node(nodes2, faces2, surface_temperature(shape2))
    end
    
    scene = poly(nodes1, faces1, color=color1, colormap=:vik, strokecolor=:black, strokewidth=1)
    poly!(nodes2, faces2, color=color2, colormap=:vik, strokecolor=:black, strokewidth=1)
    display(scene)
end



"""
"""
# function draw(shape1::ShapeModel, shape2::ShapeModel, savepath)

#     faces = VectorVector2Matrix(shape1.faces)
#     nodes = VectorVector2Matrix(shape1.nodes)

#     # limits = (-2000, 2000, -2000, 2000, -2000, 2000)
#     # aspect = (1, 1, 1)
#     fig, ax, l = poly(nodes, faces,
#         color=:gray, strokecolor=:black, strokewidth=0.1, transparency=false,
#         axis = (; type=Axis3, protrusions=(0,0,0,0), viewmode=:fit, aspect=:data)
#     )

#     # fig, ax, l = lines(points,
#     #     color=colors, colormap=:inferno, transparency=true,
#     #     axis = (; type=Axis3, protrusions=(0, 0, 0, 0), viewmode=:fit, limits=(-30, 30, -30, 30, 0, 50))
#     # )

#     set_theme!(theme_black())
    
#     record(fig, savepath, 1:120) do frame
#         # for i in 1:50
#         #     push!(points[], step!(attractor))
#         #     push!(colors[], frame)
#         # end
#         # ax.azimuth[] = 1.7pi + 0.3 * sin(2pi * frame / 120)
#         ax.elevation[] = 0.3 * sin(2π * frame / 120)
#         # notify.((points, colors))
#         l.colorrange = (0, frame)
#     end
# end
