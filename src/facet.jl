

# ################################################################
# #                      Face properties
# ################################################################

face_center(vs::StaticVector{3, <:StaticVector{3}}) = face_center(vs...)
face_center(v1::StaticVector{3}, v2::StaticVector{3}, v3::StaticVector{3}) = (v1 + v2 + v3) / 3

face_normal(vs::StaticVector{3, <:StaticVector{3}}) = face_normal(vs...)
face_normal(v1::StaticVector{3}, v2::StaticVector{3}, v3::StaticVector{3}) = normalize((v2 - v1) × (v3 - v2))

face_area(vs::StaticVector{3, <:StaticVector{3}}) = face_area(vs...)
face_area(v1::StaticVector{3}, v2::StaticVector{3}, v3::StaticVector{3}) = norm((v2 - v1) × (v3 - v2)) / 2


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
#                           Raycast
################################################################

"""
    raycast(A, B, C, R) -> Bool

Intersection detection between ray R and triangle ABC.
Note that the starting point of the ray is the origin (0, 0, 0).
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

Intersection detection between ray R and triangle ABC.
Use when the starting point of the ray is an arbitrary point `O`.
"""
raycast(A::StaticVector{3}, B::StaticVector{3}, C::StaticVector{3}, R::StaticVector{3}, O::StaticVector{3}) = raycast(A - O, B - O, C - O, R)


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


