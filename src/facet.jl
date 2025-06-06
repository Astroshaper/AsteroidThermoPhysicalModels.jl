

# ################################################################
# #                      Face properties
# ################################################################

"""
    face_center(vs::StaticVector{3, <:StaticVector{3}}) -> center
    face_center(v1::StaticVector{3}, v2::StaticVector{3}, v3::StaticVector{3}) -> center

Calculate the centroid (center) of a triangular face.

# Arguments
- `vs` : A static vector containing three vertices of the triangle
- `v1`, `v2`, `v3` : Three vertices of the triangle

# Returns
- `center::StaticVector{3}` : The centroid position vector of the triangle

# Example
```julia
v1 = SVector(0.0, 0.0, 0.0)
v2 = SVector(1.0, 0.0, 0.0)
v3 = SVector(0.0, 1.0, 0.0)
center = face_center(v1, v2, v3)  # Returns SVector(1/3, 1/3, 0.0)
```
"""
face_center(vs::StaticVector{3, <:StaticVector{3}}) = face_center(vs...)
face_center(v1::StaticVector{3}, v2::StaticVector{3}, v3::StaticVector{3}) = (v1 + v2 + v3) / 3

"""
    face_normal(vs::StaticVector{3, <:StaticVector{3}}) -> n̂
    face_normal(v1::StaticVector{3}, v2::StaticVector{3}, v3::StaticVector{3}) -> n̂

Calculate the unit normal vector of a triangular face using the right-hand rule.
The normal direction follows the ordering of vertices: v1 → v2 → v3.

# Arguments
- `vs` : A static vector containing three vertices of the triangle
- `v1`, `v2`, `v3` : Three vertices of the triangle in counter-clockwise order

# Returns
- `n̂::StaticVector{3}` : The unit normal vector of the triangle

# Mathematical Formula
The normal vector is calculated as: n̂ = normalize((v2 - v1) × (v3 - v2))

# Example
```julia
v1 = SVector(0.0, 0.0, 0.0)
v2 = SVector(1.0, 0.0, 0.0)
v3 = SVector(0.0, 1.0, 0.0)
n̂ = face_normal(v1, v2, v3)  # Returns SVector(0.0, 0.0, 1.0)
```
"""
face_normal(vs::StaticVector{3, <:StaticVector{3}}) = face_normal(vs...)
face_normal(v1::StaticVector{3}, v2::StaticVector{3}, v3::StaticVector{3}) = normalize((v2 - v1) × (v3 - v2))

"""
    face_area(vs::StaticVector{3, <:StaticVector{3}}) -> area
    face_area(v1::StaticVector{3}, v2::StaticVector{3}, v3::StaticVector{3}) -> area

Calculate the area of a triangular face using the cross product method.

# Arguments
- `vs` : A static vector containing three vertices of the triangle
- `v1`, `v2`, `v3` : Three vertices of the triangle

# Returns
- `area::Float64` : The area of the triangle [unit²]

# Mathematical Formula
The area is calculated as: A = ||(v2 - v1) × (v3 - v2)|| / 2

# Example
```julia
v1 = SVector(0.0, 0.0, 0.0)
v2 = SVector(1.0, 0.0, 0.0)
v3 = SVector(0.0, 1.0, 0.0)
area = face_area(v1, v2, v3)  # Returns 0.5
```
"""
face_area(vs::StaticVector{3, <:StaticVector{3}}) = face_area(vs...)
face_area(v1::StaticVector{3}, v2::StaticVector{3}, v3::StaticVector{3}) = norm((v2 - v1) × (v3 - v2)) / 2


################################################################
#                 Face-to-face interactions
################################################################

"""
    view_factor(cᵢ, cⱼ, n̂ᵢ, n̂ⱼ, aⱼ) -> fᵢⱼ, dᵢⱼ, d̂ᵢⱼ

Calculate the view factor from facet i to facet j, assuming Lambertian emission.
The view factor represents the fraction of radiation leaving surface i that directly reaches surface j.

# Arguments
- `cᵢ::StaticVector{3}` : Center position of facet i
- `cⱼ::StaticVector{3}` : Center position of facet j
- `n̂ᵢ::StaticVector{3}` : Unit normal vector of facet i
- `n̂ⱼ::StaticVector{3}` : Unit normal vector of facet j
- `aⱼ::Float64` : Area of facet j [m²]

# Returns
- `fᵢⱼ::Float64` : View factor from facet i to facet j [-]
- `dᵢⱼ::Float64` : Distance between facet centers [m]
- `d̂ᵢⱼ::StaticVector{3}` : Unit direction vector from facet i to facet j

# Mathematical Formula
The view factor for differential areas is calculated as:
```
fᵢⱼ = (cosθᵢ × cosθⱼ) / (π × dᵢⱼ²) × aⱼ
```
where:
- θᵢ is the angle between n̂ᵢ and the line connecting the centers
- θⱼ is the angle between n̂ⱼ and the line connecting the centers

# Diagram
```
(i)   fᵢⱼ   (j)
 △    -->    △
 cᵢ          cⱼ  : Center of each face
 n̂ᵢ          n̂ⱼ  : Normal vector of each face
```

# References
- Howell, J. R., et al. (2010). Thermal Radiation Heat Transfer
"""
function view_factor(cᵢ, cⱼ, n̂ᵢ, n̂ⱼ, aⱼ)
    dᵢⱼ = norm(cⱼ - cᵢ)       # Distance from i to j
    d̂ᵢⱼ = normalize(cⱼ - cᵢ)  # Direction from i to j

    cosθᵢ = n̂ᵢ ⋅ d̂ᵢⱼ
    cosθⱼ = n̂ⱼ ⋅ (-d̂ᵢⱼ)

    fᵢⱼ = cosθᵢ * cosθⱼ / (π * dᵢⱼ^2) * aⱼ
    fᵢⱼ, dᵢⱼ, d̂ᵢⱼ
end


################################################################
#                           Raycast
################################################################

"""
    raycast(A, B, C, R) -> Bool

Detect intersection between a ray and a triangle using the Möller-Trumbore algorithm.
The ray starts from the origin (0, 0, 0) and extends in direction R.

# Arguments
- `A::StaticVector{3}` : First vertex of the triangle
- `B::StaticVector{3}` : Second vertex of the triangle
- `C::StaticVector{3}` : Third vertex of the triangle
- `R::StaticVector{3}` : Direction vector of the ray (does not need to be normalized)

# Returns
- `Bool` : `true` if the ray intersects the triangle, `false` otherwise

# Algorithm
Uses the Möller-Trumbore ray-triangle intersection algorithm which:
1. Computes the intersection point using barycentric coordinates
2. Checks if the intersection point lies within the triangle
3. Ensures the intersection occurs in the positive ray direction (t > 0)

# Example
```julia
A = SVector(1.0, 0.0, 0.0)
B = SVector(0.0, 1.0, 0.0)
C = SVector(0.0, 0.0, 1.0)
R = SVector(1.0, 1.0, 1.0)  # Ray direction
intersects = raycast(A, B, C, R)
```

# References
- Möller, T., & Trumbore, B. (1997). Fast, minimum storage ray-triangle intersection
"""
function raycast(A::StaticVector{3}, B::StaticVector{3}, C::StaticVector{3}, R::StaticVector{3})
    E1 = B - A
    E2 = C - A
    T  = - A

    P = R × E2
    Q = T × E1
    
    P_dot_E1 = P ⋅ E1
        
    u = (P ⋅ T)  / P_dot_E1
    v = (Q ⋅ R)  / P_dot_E1
    t = (Q ⋅ E2) / P_dot_E1

    return 0 ≤ u ≤ 1 && 0 ≤ v ≤ 1 && 0 ≤ u + v ≤ 1 && t > 0
end

"""
    raycast(A, B, C, R, O) -> Bool

Detect intersection between a ray and a triangle, with arbitrary ray origin.
This is a convenience wrapper that translates the triangle to place the ray origin at (0,0,0).

# Arguments
- `A::StaticVector{3}` : First vertex of the triangle
- `B::StaticVector{3}` : Second vertex of the triangle
- `C::StaticVector{3}` : Third vertex of the triangle
- `R::StaticVector{3}` : Direction vector of the ray (does not need to be normalized)
- `O::StaticVector{3}` : Origin point of the ray

# Returns
- `Bool` : `true` if the ray intersects the triangle, `false` otherwise

# Example
```julia
A = SVector(1.0, 0.0, 0.0)
B = SVector(0.0, 1.0, 0.0)
C = SVector(0.0, 0.0, 1.0)
O = SVector(0.5, 0.5, -1.0)  # Ray origin
R = SVector(0.0, 0.0, 1.0)   # Ray direction
intersects = raycast(A, B, C, R, O)
```
"""
raycast(A::StaticVector{3}, B::StaticVector{3}, C::StaticVector{3}, R::StaticVector{3}, O::StaticVector{3}) = raycast(A - O, B - O, C - O, R)


"""
    find_visiblefacets!(shape::ShapeModel; show_progress=true)

Find facets that is visible from the facet where the observer is located.

# Arguments
- `shape` : Shape model of an asteroid

# Keyword arguments
- `show_progress` : Switch to show a progress meter
"""
function find_visiblefacets!(shape::ShapeModel; show_progress=true)
    nodes = shape.nodes
    faces = shape.faces
    face_centers = shape.face_centers
    face_normals = shape.face_normals
    face_areas = shape.face_areas
    visiblefacets = shape.visiblefacets

    ## `ProgressMeter` setting
    if show_progress
        p = Progress(length(faces); dt=1, desc="Searching for visible faces...", showspeed=true)
        ProgressMeter.ijulia_behavior(:clear)
    end

    for i in eachindex(faces)
        cᵢ = face_centers[i]
        n̂ᵢ = face_normals[i]
        aᵢ = face_areas[i]

        candidates = Int64[]
        for j in eachindex(faces)
            i == j && continue
            cⱼ = face_centers[j]
            n̂ⱼ = face_normals[j]

            Rᵢⱼ = cⱼ - cᵢ                                         # Vector from facet i to j
            Rᵢⱼ ⋅ n̂ᵢ > 0 && Rᵢⱼ ⋅ n̂ⱼ < 0 && push!(candidates, j)  # if two faces are facing each other
        end
        
        for j in candidates
            j in (visiblefacet.id for visiblefacet in visiblefacets[i]) && continue
            cⱼ = face_centers[j]
            n̂ⱼ = face_normals[j]
            aⱼ = face_areas[j]

            Rᵢⱼ = cⱼ - cᵢ    # Vector from facet i to j
            dᵢⱼ = norm(Rᵢⱼ)  # Distance from facet i to j
            
            blocked = false
            for k in candidates
                j == k && continue
                cₖ = face_centers[k]

                Rᵢₖ = cₖ - cᵢ    # Vector from facet i to k
                dᵢₖ = norm(Rᵢₖ)  # Distance from facet i to k
                
                dᵢⱼ < dᵢₖ && continue
                
                if raycast(nodes[faces[k]]..., Rᵢⱼ, cᵢ)  # if facet k blocks the view to facet j
                    blocked = true
                    break
                end
            end

            blocked && continue
            push!(visiblefacets[i], VisibleFacet(j, view_factor(cᵢ, cⱼ, n̂ᵢ, n̂ⱼ, aⱼ)...))  # i -> j
            push!(visiblefacets[j], VisibleFacet(i, view_factor(cⱼ, cᵢ, n̂ⱼ, n̂ᵢ, aᵢ)...))  # j -> i
        end

        ## Update the progress meter
        if show_progress
            showvalues = [
                ("Face ID       ", i),
                ("Visible faces ", length(visiblefacets[i])),
            ]
            ProgressMeter.next!(p; showvalues)
        end
    end
end


"""
    isilluminated(shape::ShapeModel, r☉::StaticVector{3}, i::Integer) -> Bool

Return if the `i`-th face of the `shape` model is illuminated by the direct sunlight or not

# Arguments
- `shape` : Shape model of an asteroid
- `r☉`    : Sun's position in the asteroid-fixed frame, which doesn't have to be normalized.
- `i`     : Index of the face to be checked
"""
function isilluminated(shape::ShapeModel, r☉::StaticVector{3}, i::Integer)
    cᵢ = shape.face_centers[i]
    n̂ᵢ = shape.face_normals[i]
    r̂☉ = normalize(r☉)

    n̂ᵢ ⋅ r̂☉ < 0 && return false

    for visiblefacet in shape.visiblefacets[i]
        A, B, C = shape.nodes[shape.faces[visiblefacet.id]]
        raycast(A, B, C, r̂☉, cᵢ) && return false
    end
    return true
end


