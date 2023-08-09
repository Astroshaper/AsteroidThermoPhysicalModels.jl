

"""
    grid_to_facets(xs::AbstractVector, ys::AbstractVector, zs::AbstractMatrix) -> nodes, faces

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
function grid_to_facets(xs::AbstractVector, ys::AbstractVector, zs::AbstractMatrix)
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
            
            push!(faces, ABC)
            push!(faces, DCB)
        end
    end

    return nodes, faces
end


# ################################################################
# #                      Face properties
# ################################################################

face_center(vs) = face_center(vs...)
face_center(v1, v2, v3) = (v1 + v2 + v3) / 3

face_normal(vs) = face_normal(vs...)
face_normal(v1, v2, v3) = normalize((v2 - v1) × (v3 - v2)) 

face_area(vs) = face_area(vs...)
face_area(v1, v2, v3) = norm((v2 - v1) × (v3 - v2)) / 2


################################################################
#                 Face-to-face interactions
################################################################

"""
    view_factor(cᵢ, cⱼ, n̂ᵢ, n̂ⱼ, aⱼ) -> fᵢⱼ, dᵢⱼ, d̂ᵢⱼ

View factor from facet i to j, assuming Lambertian emission.

- ---------------
- (i)   fᵢⱼ   (j)
-  △    -->    △
- ---------------
-  cᵢ          cⱼ  : Center of each face
-  n̂ᵢ          n̂ⱼ  : Normal vector of each face
-  -           aⱼ  : Area of j-th face
- ---------------
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
#                           Orinet3D
################################################################

"""
    isAbove(A, B, C, D) -> Bool

Determine if point D is above triangle face ABC.
"""
function isAbove(A, B, C, D)
    G = SA_F64[
        A[1]-D[1] A[2]-D[2] A[3]-D[3]
        B[1]-D[1] B[2]-D[2] B[3]-D[3]
        C[1]-D[1] C[2]-D[2] C[3]-D[3]
    ]

    return det(G) < 0
end

"""
    isBelow(A, B, C, D) -> Bool

Determine if point D is below triangle face ABC.
"""
function isBelow(A, B, C, D)
    G = SA_F64[
        A[1]-D[1] A[2]-D[2] A[3]-D[3]
        B[1]-D[1] B[2]-D[2] B[3]-D[3]
        C[1]-D[1] C[2]-D[2] C[3]-D[3]
    ]

    return det(G) > 0
end


################################################################
#                           Raycast
################################################################

"""
    raycast(A, B, C, R) -> Bool

Intersection detection between ray R and triangle ABC.
Note that the starting point of the ray is the origin (0, 0, 0).
"""
function raycast(A, B, C, R)
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

Intersection detection between ray R and triangle ABC.
Use when the starting point of the ray is an arbitrary point `O`.
"""
raycast(A, B, C, R, O) = raycast(A - O, B - O, C - O, R)


"""
    find_visiblefacets!(obs::Facet, facets)

Find facets that is visible from the facet where the observer is located.

# Parameters
- `obs`    : Facet where the observer stands
- `facets` : Array of `Facet`
"""
function find_visiblefacets!(shape::ShapeModel)
    nodes = shape.nodes
    faces = shape.faces
    face_centers = shape.face_centers
    face_normals = shape.face_normals
    face_areas = shape.face_areas
    visiblefacets = shape.visiblefacets

    @showprogress 1 "Searching for visible faces..." for i in eachindex(faces)
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
    end
end


"""
    isilluminated(shape::ShapeModel, r☉::AbstractVector, i::Integer) -> Bool

Return if the `i`-th face of the shape model is illuminated by the direct sunlight or not
"""
function isilluminated(shape::ShapeModel, r☉::AbstractVector, i::Integer)
    nodes = shape.nodes
    faces = shape.faces
    cᵢ = shape.face_centers[i]
    n̂ᵢ = shape.face_normals[i]
    r̂☉ = normalize(r☉)

    n̂ᵢ ⋅ r̂☉ < 0 && return false

    for visiblefacet in shape.visiblefacets[i]
        A, B, C = nodes[faces[visiblefacet.id]]
        raycast(A, B, C, r̂☉, cᵢ) && return false
    end
    return true
end


