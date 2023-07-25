"""
    mutable struct Flux

Energy flux to a facet

# Fields
- `sun`  : Flux of solar radiation,    F_sun
- `scat` : Flux of scattered sunlight, F_scat
- `rad`  : Flux of thermal radiation,  F_rad
"""
mutable struct Flux
    sun ::Float64
    scat::Float64
    rad ::Float64
end

Flux() = Flux(0., 0., 0.)

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

################################################################
#                  Triangular surface facet
################################################################

"""
    struct Facet

Triangular surface facet of a polyhedral shape model.

Note that the mesh normal indicates outward the polyhedron.

# Fields
- `center` : Position of mesh center
- `normal` : Normal vector to mesh
- `area  `   : Area of mesh
    
- `visiblefacets` : 1-D array of `VisibleFacet`
- `flux         ` : Energy flux from surrounding facets
- `temps        ` : Temperature profile in depth direction
- `_temps_      ` : Pre-allocated vector for updating temperature profile
- `force        ` : Photon recoil force
"""
struct Facet
    center::SVector{3, Float64}
    normal::SVector{3, Float64}
    area  ::Float64
    
    visiblefacets::Vector{VisibleFacet}
    flux         ::Flux
    temps        ::Vector{Float64}
    _temps_      ::Vector{Float64}
    force        ::Vector{Float64}
end

Facet(vs) = Facet(
    facet_center(vs), facet_normal(vs), facet_area(vs),
    VisibleFacet[], Flux(), Float64[], Float64[], zeros(3)
)


function Base.show(io::IO, facet::Facet)
    msg = "Surface facet\n"
    msg *= "-------------\n"

    msg *= "Center : $(facet.center)\n"
    msg *= "Normal : $(facet.normal)\n"
    msg *= "Area   : $(facet.area)\n"
    
    if isempty(facet.visiblefacets)
        msg *= "No visible facets from this facet.\n"
    else
        msg *= "$(length(facet.visiblefacets)) facets are visible from this facet:\n"
        df = DataFrame(
            id = [visiblefacet.id for visiblefacet in facet.visiblefacets],
            f  = [visiblefacet.f  for visiblefacet in facet.visiblefacets],
            d  = [visiblefacet.d  for visiblefacet in facet.visiblefacets],
        )
        msg *= "$(df)\n"
    end
    print(io, msg)
end


"""
    grid_to_facets(xs::AbstractVector, ys::AbstractVector, zs::AbstractMatrix) -> nodes, faces, facets

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
    facets = Facet[]

    for y in ys
        for x in xs
            push!(nodes, [x, y, 0])
        end
    end

    for j in eachindex(ys)[begin:end-1]
        for i in eachindex(xs)[begin:end-1]
            ABC = [i + (j-1)*length(xs), i+1 + (j-1)*length(xs), i + j*length(xs)]  # Indices of nodes of a facet ABC
            DCB = [i+1 + j*length(xs), i + j*length(xs), i+1 + (j-1)*length(xs)]    # Indices of nodes of a facet DCB
            
            push!(faces, ABC)
            push!(faces, DCB)

            A = SVector{3, Float64}(xs[i  ], ys[j  ], zs[i  , j  ])
            B = SVector{3, Float64}(xs[i+1], ys[j  ], zs[i+1, j  ])
            C = SVector{3, Float64}(xs[i  ], ys[j+1], zs[i  , j+1])
            D = SVector{3, Float64}(xs[i+1], ys[j+1], zs[i+1, j+1])

            push!(facets, Facet((A, B, C)))
            push!(facets, Facet((D, C, B)))
        end
    end
    nodes, faces, facets
end


# ################################################################
# #                      Facet properties
# ################################################################

facet_center(vs) = facet_center(vs...)
facet_center(v1, v2, v3) = (v1 + v2 + v3) / 3

facet_normal(vs) = facet_normal(vs...)
facet_normal(v1, v2, v3) = normalize((v2 - v1) × (v3 - v2)) 

facet_area(vs) = facet_area(vs...)
facet_area(v1, v2, v3) = norm((v2 - v1) × (v3 - v2)) / 2


################################################################
#                 Face-to-face interactions
################################################################


VisibleFacet(i::Facet, j::Facet, id) = VisibleFacet(id, view_factor(i, j)...)


"""
    view_factor(i::Facet, j::Facet) -> fᵢⱼ, dᵢⱼ, d̂ᵢⱼ

View factor from facet i to j, assuming Lambertian emission
"""
function view_factor(i::Facet, j::Facet)
    d⃗ᵢⱼ = j.center - i.center  # vector from facet i to j
    dᵢⱼ = norm(d⃗ᵢⱼ)
    d̂ᵢⱼ = normalize(d⃗ᵢⱼ)

    cosθᵢ = i.normal ⋅ d̂ᵢⱼ
    cosθⱼ = j.normal ⋅ (-d̂ᵢⱼ)

    fᵢⱼ = view_factor(cosθᵢ, cosθⱼ, dᵢⱼ, j.area)
    fᵢⱼ, dᵢⱼ, d̂ᵢⱼ
end

view_factor(cosθᵢ, cosθⱼ, dᵢⱼ, aⱼ) = cosθᵢ * cosθⱼ / (π * dᵢⱼ^2) * aⱼ


################################################################
#                           Orinet3D
################################################################

"""
    isFace(obs::Facet, tar::Facet) -> true/false

Determine if the two facets are facing each other
"""
isFace(obs::Facet, tar::Facet) = (tar.center - obs.center) ⋅ tar.normal < 0

"""
    isAbove(A, B, C, D)             -> Bool

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
    isBelow(A, B, C, D)             -> Bool

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
function find_visiblefacets!(nodes, faces, facets)
    for i in eachindex(faces)
        candidates = Int64[]
        for j in eachindex(faces)
            i == j && continue
            isAbove(nodes[faces[i]]..., facets[j].center) && isFace(facets[i], facets[j]) && push!(candidates, j)
        end
        
        for j in candidates
            j in (visiblefacet.id for visiblefacet in facets[i].visiblefacets) && continue

            Rⱼ = facets[j].center - facets[i].center  # Vector from facet i to j
            dⱼ = norm(Rⱼ)                             # Distance from facet i to j
            
            blocked = false
            for k in candidates
                j == k && continue
                Rₖ = facets[k].center - facets[i].center  # Vector from facet i to k
                dₖ = norm(Rₖ)                             # Distance from facet i to k
                
                dⱼ < dₖ && continue
                
                if raycast(nodes[faces[k]]..., Rⱼ, facets[i].center)      # if facet k blocks the view to facet j
                    blocked = true
                    break
                end
            end

            blocked && continue
            push!(facets[i].visiblefacets, VisibleFacet(facets[i], facets[j], j))
            push!(facets[j].visiblefacets, VisibleFacet(facets[j], facets[i], i))
        end
    end
end


# """
#     This function will be reused when parallelizing the code.
# """
# function find_visiblefacets!(obs::Facet, facets)
    
#     candidates = Int64[]
#     for (j, facet) in enumerate(facets)
#         isAbove(obs, facet) && isFace(obs, facet) && push!(candidates, j)
#     end

#     for j in candidates
#         Rⱼ = facets[j].center - obs.center   # Vector from the facet `obs` to j
#         dⱼ = norm(Rⱼ)                        # Distance to facet j

#         blocked = false
#         for k in candidates
#             j == k && continue
#             Rₖ = facets[k].center - obs.center  # Vector from the facet `obs` to k
#             dₖ = norm(Rₖ)                       # Distance to facet k

#             dⱼ < dₖ && continue

#             if raycast(facets[k], Rⱼ, obs)      # Not visible because the facet k blocks the view to j
#                 blocked = true
#                 break
#             end
#         end
#         blocked && continue

#         push!(obs.visiblefacets, VisibleFacet(obs, facets[j], j))
#     end
# end

# """
#     find_visiblefacets!(facets)

# Find facets that is visible from each facet

# # Parameters
# - `facets` : Array of `Facet`
# """
# function find_visiblefacets!(facets)
#     for obs in facets
#         find_visiblefacets!(obs, facets)
#     end
# end

"""
    isIlluminated(obs::Facet, r̂☉, facets) -> Bool

Return if the observation facet is illuminated by the direct sunlight or not
"""
function isIlluminated(obs::Facet, r̂☉, shape)
    obs.normal ⋅ r̂☉ < 0 && return false
    for visiblefacet in obs.visiblefacets
        A, B, C = shape.nodes[shape.faces[visiblefacet.id]]
        raycast(A, B, C, r̂☉, obs.center) && return false
    end
    return true
end


