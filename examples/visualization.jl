using GLMakie, CairoMakie
using ScatteredInterpolation

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
function draw(shape::ShapeModel; data=nothing, r̂☉=[1,0,0.], colormap=:viridis, strokecolor=:gray20, strokewidth=1)
    nodes = VectorVector2Matrix(shape.nodes)
    faces = VectorVector2Matrix(shape.faces)

    set_theme!(backgroundcolor=:black)
    
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
        strokecolor=strokecolor, strokewidth=strokewidth,
        size=(1500,1500)
    )

    set_theme!(backgroundcolor=:white)

    display(scene)

    # fig = Figure()
    # ax = Axis3(fig[1, 1], aspect=:data)

    # poly!(ax, nodes, faces,
    #     color=color, colormap=colormap,
    #     strokecolor=strokecolor, strokewidth=1,
    #     size=(1500,1500)
    # )

    # Colorbar(fig[1, 2], limits = (minimum(surf_temps), maximum(surf_temps)), colormap=colormap, ticks = 0:50:400)

    # display(fig)
end


function draw(shapes::Tuple, sec_from_pri, R₂₁; data=nothing, r̂☉=[1,0,0.], colormap=:viridis, strokecolor=:gray20, strokewidth=1)

    fig = Figure(backgroundcolor=:gray40)
    ax = Axis3(fig[1, 1], aspect=:data)

    for (i, shape) in enumerate(shapes)

        i == 1 && (nodes = shape.nodes)
        i == 2 && (nodes = [sec_from_pri + R₂₁ * node for node in shape.nodes])
        
        nodes = VectorVector2Matrix(nodes)
        faces = VectorVector2Matrix(shape.faces)
        
        color = face2node(nodes, faces, surface_temperature(shape))
        poly!(ax, nodes, faces; colormap, color)

        if i == 2
            limits = extrema(vcat(surface_temperature.(shapes)...))
            Colorbar(fig[1, 2]; colormap, limits, label="Temperature [K]")
        end
    end
    display(fig)
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


################################################################
#                      Temperature map
################################################################


latitude(facet::Facet) = latitude(facet.center)
latitude(r) = asin(r[3] / norm(r))

longitude(facet::Facet) = longitude(facet.center)
longitude(r) = atan(r[2], r[1])

"""
    facet_to_grid(shape, data) -> x, y, gridded

Make a lat-lon grid from facet-correlated data using ScatteredInterpolation.jl
"""
function facet_to_grid(shape, data)

    lons = rad2deg.(longitude(facet) for facet in shape.facets)
    lats = rad2deg.(latitude(facet)  for facet in shape.facets)
    points = hcat(lons, lats)'

    n = 180
    x = range(-180, 180, length=n)
    y = range(-90, 90, length=n)
    X = repeat(x, n)[:]
    Y = repeat(y', n)[:]
    gridPoints = [X Y]'

    itp = interpolate(Multiquadratic(), points, data)
    interpolated = evaluate(itp, gridPoints)
    gridded = reshape(interpolated, n, n)

    x, y, gridded
end


"""
    temperature_map(shape, temps=surface_temperature(shape))

Make a global 2D-map from temperature based on every facet.

ScatteredInterpolation.jl works well.
- https://eljungsk.github.io/ScatteredInterpolation.jl/dev/

CairoMakie.tricontourf seems unavailable now?
- https://docs.makie.org/v0.18.0/examples/plotting_functions/tricontourf/
"""
function temperature_map(shape::ShapeModel;
    temps=surface_temperature(shape),
    colormap=:hot, colorrange=extrema(temps),
    draw_contour=true, nlevels=15, ticks=0:20:5000,
    filepath="temp_map.png", title="",
)

    T_min, T_max = extrema(temps)
    println("Max. temperature: ", T_max)
    println("Min. temperature: ", T_min)
    
    fig = Figure(resolution=(800,500))
    ax = Axis(fig[1, 1],
        title=title,
        xlabel="Longitude [deg]",
        ylabel="Latitude [deg]",
        xticks=-180:30:180,
        yticks=-90:30:90,
    )
    xlims!(ax, -180, 180)
    ylims!(ax, -90, 90)

    x, y, gridded = facet_to_grid(shape, temps)
    levels = range(colorrange..., nlevels+1)
    cntrf = contourf!(ax, x, y, gridded; colormap, levels, extendlow=:auto, extendhigh=:auto)
    draw_contour && contour!(x, y, gridded; color=:black, linewidth=0.5, levels)
    Colorbar(fig[:, end+1], cntrf, ticks=ticks, label="Temperature [K]")

    save(filepath, fig)
    fig
end

"""
    temperature_map(shape, temps=surface_temperature(shape))

Make temperature maps of a binary asteroid, `shape1` and `shape2`.
"""
function temperature_map(shape1::ShapeModel, shape2::ShapeModel;
    temps1=surface_temperature(shape1), temps2=surface_temperature(shape2),
    colormap=:hot, colorrange=extrema(vcat(temps1, temps2)),
    draw_contour=true, nlevels=15, ticks=0:20:5000,
    filepath="temp_maps.pdf", titles=("", ""),
)

    T_min, T_max = extrema(vcat(temps1, temps2))
    println("Max. temperature: ", T_max)
    println("Min. temperature: ", T_min)

    fig = Figure(resolution=(1400, 600))

    for (idx_shape, (shape, temps, title)) in enumerate(zip((shape1, shape2), (temps1, temps2), titles))

        ax = Axis(fig[1, idx_shape],
            title=title,
            xlabel="Longitude [deg]",
            ylabel="Latitude [deg]",
            xticks=-180:30:180,
            yticks=-90:30:90,
        )
        xlims!(ax, -180, 180)
        ylims!(ax, -90, 90)
    
        x, y, gridded = facet_to_grid(shape, temps)
        levels = range(colorrange..., nlevels+1)
        cntrf = contourf!(ax, x, y, gridded; colormap, levels, extendlow=:auto, extendhigh=:auto)
        draw_contour && contour!(x, y, gridded; color=:black, linewidth=0.5, levels)

        if idx_shape == 2
            Colorbar(fig[2, :], cntrf, ticks=ticks, label="Temperature [K]", vertical=false, flipaxis=false)
        end
    end

    save(filepath, fig)
    fig
end
