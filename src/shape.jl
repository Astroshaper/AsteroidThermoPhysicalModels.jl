
"""
    struct VisibleFacet

Index of an interfacing facet and its view factor

# Fields
- `id` : Index of the interfacing facet
- `f`  : View factor from facet i to j
- `d`  : Distance from facet i to j
- `d̂`  : Normal vector from facet i to j
"""
struct VisibleFacet
    id::Int64
    f ::Float64
    d ::Float64
    d̂ ::SVector{3, Float64}
end


"""
    ShapeModel

A polyhedral shape model of an asteroid.

# Fields
- `nodes`         : Vector of node positions
- `faces`         : Vector of vertex indices of faces
- `face_centers`  : Center position of each face
- `face_normals`  : Normal vector of each face
- `face_areas`    : Area of of each face
- `visiblefacets` : Vector of vector of `VisibleFacet`
"""
mutable struct ShapeModel
    nodes        ::Vector{SVector{3, Float64}}
    faces        ::Vector{SVector{3, Int}}

    face_centers ::Vector{SVector{3, Float64}}
    face_normals ::Vector{SVector{3, Float64}}
    face_areas   ::Vector{Float64}

    visiblefacets::Vector{Vector{VisibleFacet}}
end


function Base.show(io::IO, shape::ShapeModel)
    msg = "Shape model\n"
    msg *= "-----------\n"
    msg *= "Number of nodes   : $(length(shape.nodes))\n"
    msg *= "Number of faces   : $(length(shape.faces))\n"
    msg *= "Volume            : $(polyhedron_volume(shape))\n"
    msg *= "Equivalent radius : $(equivalent_radius(shape))\n"
    msg *= "Maximum radius    : $(maximum_radius(shape))\n"
    msg *= "Minimum radius    : $(minimum_radius(shape))\n"
    print(io, msg)
end

"""
    load_shape_obj(shapepath; scale=1.0, find_visible_facets=false; show_progress=true)

Load a shape model from a Wavefront OBJ file.

# Arguments
- `shapepath` : Path to a Wavefront OBJ file

# Keyword arguments
- `scale`               : Scale factor of the shape model
- `find_visible_facets` : Switch to find visible facets
- `show_progress`       : Switch to show a progress meter
"""
function load_shape_obj(shapepath; scale=1.0, find_visible_facets=false, show_progress=true)
    # TODO: use MeshIO.jl
    nodes, faces = loadobj(shapepath; scale = scale, message = false)

    face_centers = [face_center(nodes[face]) for face in faces]
    face_normals = [face_normal(nodes[face]) for face in faces]
    face_areas   = [face_area(nodes[face])   for face in faces]

    visiblefacets = [VisibleFacet[] for _ in faces]

    shape = ShapeModel(nodes, faces, face_centers, face_normals, face_areas, visiblefacets)
    find_visible_facets && find_visiblefacets!(shape; show_progress)
    
    return shape
end


################################################################
#               Create a shape model from grid
################################################################


"""
    grid_to_faces(xs::AbstractVector, ys::AbstractVector, zs::AbstractMatrix) -> nodes, faces

Convert a regular grid (x, y) and corresponding z-coordinates to triangular facets

    | ⧹| ⧹| ⧹|
j+1 ・--C--D--・
    |⧹ |⧹ |⧹ |
    | ⧹| ⧹| ⧹|
j   ・--A--B--・
    |⧹ |⧹ |⧹ |
       i  i+1

# Arguments
- `xs::AbstractVector` : x-coordinates of grid points (should be sorted)
- `ys::AbstractVector` : y-coordinates of grid points (should be sorted)
- `zs::AbstractMatrix` : z-coordinates of grid points
"""
function grid_to_faces(xs::AbstractVector, ys::AbstractVector, zs::AbstractMatrix)
    nodes = SVector{3, Float64}[]
    faces = SVector{3, Int}[]

    for j in eachindex(ys)
        for i in eachindex(xs)
            push!(nodes, @SVector [xs[i], ys[j], zs[i, j]])
        end
    end

    for j in eachindex(ys)[begin:end-1]
        for i in eachindex(xs)[begin:end-1]
            ABC = @SVector [i + (j-1)*length(xs), i+1 + (j-1)*length(xs), i + j*length(xs)]  # Indices of nodes of △ABC
            DCB = @SVector [i+1 + j*length(xs), i + j*length(xs), i+1 + (j-1)*length(xs)]    # Indices of nodes of △DCB

            push!(faces, ABC, DCB)
        end
    end

    return nodes, faces
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
    nodes, faces = grid_to_faces(xs, ys, zs)
    nodes .*= scale
    
    face_centers = [face_center(nodes[face]) for face in faces]
    face_normals = [face_normal(nodes[face]) for face in faces]
    face_areas   = [face_area(nodes[face])   for face in faces]

    visiblefacets = [VisibleFacet[] for _ in faces]

    shape = ShapeModel(nodes, faces, face_centers, face_normals, face_areas, visiblefacets)
    find_visible_facets && find_visiblefacets!(shape)
    
    return shape
end


################################################################
#                      Shape properites
################################################################

"""
    polyhedron_volume(nodes, faces) -> vol
    polyhedron_volume(shape::ShapeModel) -> vol

Calculate the volume of a closed polyhedral shape using the divergence theorem.

# Arguments
- `nodes::Vector{SVector{3}}` : Vector of vertex positions
- `faces::Vector{SVector{3,Int}}` : Vector of triangular face definitions (vertex indices)
- `shape::ShapeModel` : Complete shape model containing nodes and faces

# Returns
- `vol::Float64` : Volume of the polyhedron [unit³]

# Algorithm
Uses the divergence theorem to compute volume from surface triangles:
```
V = (1/6) Σᵢ (Aᵢ × Bᵢ) ⋅ Cᵢ
```
where Aᵢ, Bᵢ, Cᵢ are the three vertices of triangle i.

# Requirements
- The polyhedron must be closed (watertight)
- Face normals should point outward for positive volume
- Faces must be consistently oriented (all clockwise or all counter-clockwise)

# Example
```julia
# Create a tetrahedron
nodes = [SVector(0.0, 0.0, 0.0),
         SVector(1.0, 0.0, 0.0),
         SVector(0.0, 1.0, 0.0),
         SVector(0.0, 0.0, 1.0)]
faces = [SVector(1, 2, 3),
         SVector(1, 2, 4),
         SVector(1, 3, 4),
         SVector(2, 3, 4)]
vol = polyhedron_volume(nodes, faces)  # Returns 1/6
```

# See Also
- `equivalent_radius` to convert volume to equivalent spherical radius
"""
function polyhedron_volume(nodes, faces)
    volume = 0.
    for face in faces
        A, B, C = nodes[face]
        volume += (A × B) ⋅ C / 6
    end
    volume
end

polyhedron_volume(shape::ShapeModel) = polyhedron_volume(shape.nodes, shape.faces)

"""
    equivalent_radius(VOLUME::Real) -> R_eq
    equivalent_radius(shape::ShapeModel) -> R_eq

Calculate the radius of a sphere with the same volume as the given volume or shape.

# Arguments
- `VOLUME::Real` : Volume of the object [unit³]
- `shape::ShapeModel` : Shape model to compute equivalent radius for

# Returns
- `R_eq::Float64` : Equivalent spherical radius [unit]

# Formula
```
R_eq = (3V/4π)^(1/3)
```

# Example
```julia
# For a cube with side length 2
volume = 8.0  # 2³
R_eq = equivalent_radius(volume)  # Returns ≈ 1.24

# For a shape model
shape = load_shape_obj("asteroid.obj")
R_eq = equivalent_radius(shape)
```

# See Also
- `polyhedron_volume` for volume calculation
- `maximum_radius`, `minimum_radius` for other size metrics
"""
equivalent_radius(VOLUME::Real) = (3VOLUME/4π)^(1/3)
equivalent_radius(shape::ShapeModel) = equivalent_radius(polyhedron_volume(shape))

"""
    maximum_radius(nodes) -> R_max
    maximum_radius(shape::ShapeModel) -> R_max

Find the maximum distance from the origin to any vertex in the shape.

# Arguments
- `nodes::Vector{SVector{3}}` : Vector of vertex positions
- `shape::ShapeModel` : Shape model

# Returns
- `R_max::Float64` : Maximum radius [unit]

# Notes
- Assumes the shape is centered at the origin
- Useful for determining bounding spheres
- For binary systems, ensure each component is centered appropriately

# See Also
- `minimum_radius` for the minimum distance
- `equivalent_radius` for volume-based radius
"""
maximum_radius(nodes) = maximum(norm, nodes)
maximum_radius(shape::ShapeModel) = maximum_radius(shape.nodes)

"""
    minimum_radius(nodes) -> R_min
    minimum_radius(shape::ShapeModel) -> R_min

Find the minimum distance from the origin to any vertex in the shape.

# Arguments
- `nodes::Vector{SVector{3}}` : Vector of vertex positions
- `shape::ShapeModel` : Shape model

# Returns
- `R_min::Float64` : Minimum radius [unit]

# Notes
- Assumes the shape is centered at the origin
- For highly irregular shapes, this may be much smaller than the equivalent radius
- Useful for determining the deepest surface features

# See Also
- `maximum_radius` for the maximum distance
- `equivalent_radius` for volume-based radius
"""
minimum_radius(nodes) = minimum(norm, nodes)
minimum_radius(shape::ShapeModel) = minimum_radius(shape.nodes)

