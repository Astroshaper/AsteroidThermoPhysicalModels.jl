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
- `A::T1` : Position of 1st vertex
- `B::T1` : Position of 2nd vertex
- `C::T1` : Position of 3rd vertex

- `center::T1` : Position of mesh center
- `normal::T1` : Normal vector to mesh
- `area  ::T2`   : Area of mesh
    
- `visiblefacets::T3` : 1-D array of `VisibleFacet`
- `flux         ::T4` : Energy flux from surrounding facets
- `temps        ::T5` : Temperature profile in depth direction
- `_temps_      ::T5` : Pre-allocated vector for updating temperature profile
- `force        ::T6` : Photon recoil force
"""
struct Facet
    A::SVector{3, Float64}
    B::SVector{3, Float64}
    C::SVector{3, Float64}
    
    center::SVector{3, Float64}
    normal::SVector{3, Float64}
    area  ::Float64
    
    visiblefacets::StructVector{VisibleFacet, NamedTuple{(:id, :f, :d, :d̂), Tuple{Vector{Int64}, Vector{Float64}, Vector{Float64}, Vector{SVector{3, Float64}}}}, Int64}
    flux         ::Flux
    temps        ::Vector{Float64}
    _temps_      ::Vector{Float64}
    force        ::Vector{Float64}
end

Facet(A, B, C) = Facet([A, B, C])
Facet(vs) = Facet(
    vs[1], vs[2], vs[3],
    facet_center(vs), facet_normal(vs), facet_area(vs),
    StructArray(VisibleFacet[]), Flux(), Float64[], Float64[], zeros(3)
)

function Base.show(io::IO, facet::Facet)
    msg = "Surface facet"
    msg *= "Surface facet\n"
    msg *= "-------------\n"

    msg *= "Vertices\n"
    msg *= "    A : $(facet.A)\n"
    msg *= "    B : $(facet.B)\n"
    msg *= "    C : $(facet.C)\n"

    msg *= "Center : $(facet.center)\n"
    msg *= "Normal : $(facet.normal)\n"
    msg *= "Area   : $(facet.area)\n"
    
    if isempty(facet.visiblefacets)
        msg *= "No visible facets.\n"
    else
        @unpack visiblefacets = facet
        msg *= "$(length(visiblefacets)) facets are visible:\n"
        df = DataFrame(id=visiblefacets.id, f=visiblefacets.f, d=visiblefacets.d)
        msg *= "$(df)\n"
    end
    print(io, msg)
end

"""
Array of `Facet`, converted from arrays of nodes and faces of a shape model
"""
getfacets(nodes, faces) = StructArray([Facet(nodes[face]) for face in faces])


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
isFace(obs::Facet, tar::Facet) = (tar.center - obs.center) ⋅ tar.normal < 0 ? true : false

"""
    isAbove(A, B, C, D)             -> Bool
    isAbove(facet::Facet, D)        -> Bool
    isAbove(obs::Facet, tar::Facet) -> Bool

Determine if point D is above triangle facet ABC.
"""
function isAbove(A, B, C, D)
    G = SA_F64[
        A[1]-D[1] A[2]-D[2] A[3]-D[3]
        B[1]-D[1] B[2]-D[2] B[3]-D[3]
        C[1]-D[1] C[2]-D[2] C[3]-D[3]
    ]

    det(G) < 0 ? true : false
end

"""
    isBelow(A, B, C, D)             -> Bool
    isBelow(facet::Facet, D)        -> Bool
    isBelow(obs::Facet, tar::Facet) -> Bool

Determine if point D is below triangle facet ABC.
"""
function isBelow(A, B, C, D)
    G = SA_F64[
        A[1]-D[1] A[2]-D[2] A[3]-D[3]
        B[1]-D[1] B[2]-D[2] B[3]-D[3]
        C[1]-D[1] C[2]-D[2] C[3]-D[3]
    ]

    det(G) > 0 ? true : false
end

isAbove(facet::Facet, D) = isAbove(facet.A, facet.B, facet.C, D)
isBelow(facet::Facet, D) = isBelow(facet.A, facet.B, facet.C, D)

isAbove(obs::Facet, tar::Facet) = isAbove(obs.A, obs.B, obs.C, tar.center)
isBelow(obs::Facet, tar::Facet) = isBelow(obs.A, obs.B, obs.C, tar.center)


################################################################
#                           Raycast
################################################################

"""
    raycast(A, B, C, R)                           -> Bool
    raycast(facet::Facet, R)                      -> Bool
    raycast(facet::Facet, R, obs::AbstractVector) -> Bool
    raycast(facet::Facet, R, obs::Facet)          -> Bool

Ray-triangle intersection detection.
The ray R starts from the origin.
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

    0 ≤ u ≤ 1 && 0 ≤ v ≤ 1 && 0 ≤ u + v ≤ 1 && t > 0 ? true : false
end

raycast(facet::Facet, R) = raycast(facet.A, facet.B, facet.C, R)
raycast(facet::Facet, R, obs::AbstractVector) = raycast(facet.A - obs, facet.B - obs, facet.C - obs, R)
raycast(facet::Facet, R, obs::Facet) = raycast(facet, R, obs.center)

raycast(A, B, C, R, obs::AbstractVector) = raycast(A - obs, B - obs, C - obs, R)
raycast(A, B, C, R, obs::Facet) = raycast(A, B, C, R, obs.center)


"""
    findVisibleFacets!(obs::Facet, facets)

Find facets that is visible from the facet where the observer is located.

# Parameters
- `obs`    : Facet where the observer stands
- `facets` : Array of `Facet`
"""
function findVisibleFacets!(obs::Facet, facets)
    ids = Int64[]
    for (id, tar) in enumerate(facets)
        isAbove(obs, tar) && isFace(obs, tar) && push!(ids, id)
    end
    
    ii = copy(ids)

    for i in ii
        tar_i = facets[i]
        Rᵢ = tar_i.center - obs.center
        dᵢ = norm(Rᵢ)      # distance to facet i
        for j in ii
            i == j && continue

            tar_j = facets[j]
            Rⱼ = tar_j.center - obs.center
            dⱼ = norm(Rⱼ)  # distance to facet j

            raycast(tar_j, Rᵢ, obs) && (dᵢ < dⱼ ? filter!(x->x≠j, ids) : filter!(x->x≠i, ids))
        end
    end
    
    for id in ids
        tar = facets[id]
        push!(obs.visiblefacets, VisibleFacet(obs, tar, id))
    end
end

"""
    findVisibleFacets!(facets)

Find facets that is visible from each facet

# Parameters
- `facets` : Array of `Facet`
"""
function findVisibleFacets!(facets)
    for obs in facets
        findVisibleFacets!(obs, facets)
    end
end

"""
    isIlluminated(obs::AbstractVector, r̂☉, facets) -> Bool
    isIlluminated(obs::Facet, r̂☉, facets)          -> Bool

Return if the observation point/facet is illuminated by the direct sunlight or not
"""
function isIlluminated(obs::AbstractVector, r̂☉, facets)
    for facet in facets
        raycast(facet, r̂☉, obs) && return false
    end
    return true
end

function isIlluminated(obs::Facet, r̂☉, facets)
    obs.normal ⋅ r̂☉ < 0 && return false
    for id in obs.visiblefacets.id
        raycast(facets[id], r̂☉, obs) && return false
    end
    return true
end

"""
    isAboveHorizon(facet::Facet) -> Bool

Return if the facet is above its local horizon or not
"""
isAboveHorizon(facet::Facet) = length(facet.visiblefaces) == 0

