

################################################################
#                      View factor
################################################################

"""
    ViewFactor

Index of an interfacing mesh and its view factor

# Fields

- `id`  : Index of the interfacing mesh
- `fᵢⱼ` : View factor from mesh i to mesh j
"""
struct ViewFactor
    id::Int64
    fᵢⱼ::Float64
end


"""
    getViewFactor(mᵢ, mⱼ) -> fᵢⱼ

View factor from mesh i to mesh j, assuming Lambertian emission
"""
function getViewFactor(mᵢ, mⱼ)
    d⃗ᵢⱼ = mⱼ.center - mᵢ.center  # vector from mesh i to mesh j
    d̂ᵢⱼ = normalize(d⃗ᵢⱼ)
    dᵢⱼ = norm(d⃗ᵢⱼ)

    cosθᵢ = mᵢ.normal ⋅ d̂ᵢⱼ
    cosθⱼ = mⱼ.normal ⋅ (-d̂ᵢⱼ)

    fᵢⱼ = getViewFactor(cosθᵢ, cosθⱼ, dᵢⱼ, mⱼ.area)
end

getViewFactor(cosθᵢ, cosθⱼ, dᵢⱼ, aⱼ) = cosθᵢ * cosθⱼ / (π * dᵢⱼ^2) * aⱼ


################################################################
#                 Face-to-face interactions
################################################################

"""
    Flux{T}

Energy flux from/to a mesh

SMeshにfluxesを保持させる

NamedTuple : 作るのは早い (immutable)
Dict : 作るのは遅い。その後に計算で使う分には遅くないかも？ (mutable)
それならmutable structがいいか（使い方次第では、遅くならない？）
サイズが小さいなら、Dictより、mutable structが良さそう

フィールドの値を更新すると、allocationがかなり発生するかもしれない
遅くなりそうなら、Vector{Float64}で代用する
https://bkamins.github.io/julialang/2020/10/16/gctime.html

# Fields
- `sun`  : F_sun
- `scat` : F_scat
- `rad`  : F_rad
- `cond` : k(dT/dx)
- `ϵσT⁴` : ϵσT⁴

surface roughness infrared beamingの効果も実装する
"""
mutable struct Flux{T}
    sun::T
    scat::T
    rad::T
    cond::T
    ϵσT⁴::T
end


################################################################
#                  Triangular surface mesh
################################################################


"""
    SMesh{T1, T2, T3}

Triangular surface mesh of a polyhedral shape model.

Note that the mesh normal indicates outward the polyhedron.

# Fields
- `A` : Position of 1st vertex
- `B` : Position of 2nd vertex
- `C` : Position of 3rd vertex

- `center` : Position of mesh center
- `normal` : Normal vector to mesh
- `area`   : Area of mesh
    
- `viewfactors` : 1-D array of `ViewFactor`
- # fluxes::T4  # 
"""
struct SMesh{T1, T2, T3, T4}
    A::T1
    B::T1
    C::T1
    
    center::T1
    normal::T1
    area::T2
    
    viewfactors::T3
    flux::T4
end


function Base.show(io::IO, smesh::SMesh)
    println(io, "Surface mesh")
    println(io, "------------")

    println("Vertices")
    println("    A : ", smesh.A)
    println("    B : ", smesh.B)
    println("    C : ", smesh.C)

    println("Center : ", smesh.center)
    println("Normal : ", smesh.normal)
    println("Area   : ", smesh.area)

    length(smesh.viewfactors) != 0 && println("Visible faces")
    length(smesh.viewfactors) != 0 && println(smesh.viewfactors.id)
    length(smesh.viewfactors) != 0 && println(smesh.viewfactors.fᵢⱼ)
end


SMesh(A, B, C) = SMesh([A, B, C])
SMesh(vs) = SMesh(vs[1], vs[2], vs[3], getcenter(vs), getnormal(vs), getarea(vs), StructArray(ViewFactor[]), Flux(0.,0.,0.,0.,0.))

getmeshes(nodes, faces) = StructArray([SMesh(nodes[face]) for face in faces])


################################################################
#                      Mesh properites
################################################################

getcenter(vs) = getcenter(vs[1], vs[2], vs[3])
getcenter(v1, v2, v3) = (v1 + v2 + v3) / 3

getnormal(vs) = getnormal(vs[1], vs[2], vs[3])
getnormal(v1, v2, v3) = normalize((v2 - v1) × (v3 - v2)) 

getarea(vs) = getarea(vs[1], vs[2], vs[3])
getarea(v1, v2, v3) = norm((v2 - v1) × (v3 - v2)) * 0.5

getcenters(meshes) = [m.center for m in meshes]
getnormals(meshes) = [m.normal for m in meshes]
getareas(meshes) = [m.area for m in meshes]


################################################################
#                           Orinet3D
################################################################

isAbove(mesh::SMesh, D) = isAbove(mesh.A, mesh.B, mesh.C, D)
isBelow(mesh::SMesh, D) = isBelow(mesh.A, mesh.B, mesh.C, D)

isAbove(obs::SMesh, tar::SMesh) = isAbove(obs.A, obs.B, obs.C, tar.center)
isBelow(obs::SMesh, tar::SMesh) = isBelow(obs.A, obs.B, obs.C, tar.center)

function isAbove(A, B, C, D)
    G = SA_F64[
        A[1]-D[1] A[2]-D[2] A[3]-D[3]
        B[1]-D[1] B[2]-D[2] B[3]-D[3]
        C[1]-D[1] C[2]-D[2] C[3]-D[3]
    ]

    det(G) < 0 ? true : false
end


function isBelow(A, B, C, D)
    G = SA_F64[
        A[1]-D[1] A[2]-D[2] A[3]-D[3]
        B[1]-D[1] B[2]-D[2] B[3]-D[3]
        C[1]-D[1] C[2]-D[2] C[3]-D[3]
    ]

    det(G) > 0 ? true : false
end


isFace(obs::SMesh, tar::SMesh) = isFace(obs.center, tar)
isFace(obs::AbstractVector, tar::SMesh) = (tar.center - obs) ⋅ tar.normal < 0 ? true : false


################################################################
#                           Raycast
################################################################

## Viewd from an observor's mesh/position
raycast(mesh, R, obs::SMesh) = raycast(mesh, R, obs.center)
raycast(mesh, R, obs::AbstractVector) = raycast(mesh.A - obs, mesh.B - obs, mesh.C - obs, R)

## Viewd from the origin
raycast(mesh, R) = raycast(mesh.A, mesh.B, mesh.C, R)


function raycast(A, B, C, R)
    E1 = B - A
    E2 = C - A
    T  = - A

    P = R × E2
    Q = T × E1
    
    P_dot_E1 = P ⋅ E1
        
    u = (P ⋅ T)  / P_dot_E1
    v = (Q ⋅ R)  / P_dot_E1
    t = (Q ⋅ E2) / P_dot_E1  # Distance to triangle (if t < 0, the ray hits from the back)

    0 ≤ u ≤ 1 && 0 ≤ v ≤ 1 && 0 ≤ u + v ≤ 1 && t > 0 ? true : false
end


"""
    findVisibleFaces!(obs::SMesh, meshes)

Find faces directly seen from the observer on the face

# Parameters
- `obs`    : Mesh where the observer stands
- `meshes` : Array of SMesh instances
"""
function findVisibleFaces!(obs::SMesh, meshes)
    ids = Int64[]
    for i in eachindex(meshes)
        tar = meshes[i]
        isAbove(obs, tar) && isFace(obs, tar) && push!(ids, i)
    end

    for i in copy(ids)
        tar_i = meshes[i]
        R = tar_i.center - obs.center
        for j in copy(ids)
            i == j && continue
            tar_j = meshes[j]

            if raycast(tar_j, R, obs)
                dᵢ = norm(R)                          # distance to i-th mesh
                dⱼ = norm(tar_j.center - obs.center)  # distance to j-th mesh
                dᵢ < dⱼ ? filter!(x->x≠j, ids) : filter!(x->x≠i, ids)
            end
        end
    end
    
    for id in ids
        fᵢⱼ = getViewFactor(obs, meshes[id])
        push!(obs.viewfactors, ViewFactor(id, fᵢⱼ))
    end
end


"""
    findVisibleFaces!(meshes)

Find faces directly seen from the observer on each face

# Parameters
- `meshes` : Array of SMesh instances
"""
function findVisibleFaces!(meshes)
    for obs in meshes
        findVisibleFaces!(obs, meshes)
    end
end


"""
Return true if the observing mesh is illuminated by the direct sunlight, false if not
"""
function isIlluminated(obs::SMesh, r̂☉, meshes)
    for id in obs.viewfactors.id
        raycast(meshes[id], r̂☉, obs) && return false
    end
    return true
end


"""
    isAboveHorizon(mesh) -> Bool

The mesh is above its local horizon or not
"""
isAboveHorizon(mesh) = length(mesh.viewfactors) == 0



# function getIlluminatedFaces(r̂☉, meshes::Vector{SMesh})
#     illuminated = Int64[]
#     for i in eachindex(meshes)
#         obs = meshes[i]
#         Ψ = obs.normal ⋅ r̂☉  # cosine of the Sun illumination angle
#         Ψ > 0 && isIlluminated(obs, r̂☉, meshes) && push!(illuminated, i)
#     end
#     illuminated
# end


# function getIlluminatedFaces_PSEUDOCONVEX(r̂☉, meshes::Vector{SMesh})
#     illuminated = Int64[]
#     for i in eachindex(meshes)
#         obs = meshes[i]
#         Ψ = obs.normal ⋅ r̂☉  # cosine of the Sun illumination angle
#         Ψ > 0 && push!(illuminated, i)
#     end
#     illuminated
# end


################################################################
#                    Solid angle of a mesh
################################################################

"""
    getSolidAngle(mesh::SMesh, observer) -> Ω

Solid angle of a triangular mesh,
equal to area of the corresponding spherical triangle
"""
function getSolidAngle(mesh::SMesh, obs::AbstractVector)

    ## Vectors from observer to mesh vertices
    A = mesh.A - obs
    B = mesh.B - obs
    C = mesh.C - obs

    AOB = getangle(A, B)
    BOC = getangle(B, C)
    COA = getangle(C, A)
    
    Ω = getSphericalExcess(AOB, BOC, COA)
end


"""
    getSphericalExcess(a, b, c) -> E

Area of a spherical triangle
c.f. L'Huilier's Theorem

# Arguments
- `a`, `b`, `c` : side length of a spherical traiangle
"""
function getSphericalExcess(a, b, c)
    s = (a + b + c) * 0.5  # semiperimeter
        
    E = tan(s*0.5)
    E *= tan((s - a)*0.5)
    E *= tan((s - b)*0.5)
    E *= tan((s - c)*0.5)
    E = sqrt(E)
    E = 4 * atan(E)
end


getangle(v1, v2) = acos(normalize(v1) ⋅ normalize(v2))
