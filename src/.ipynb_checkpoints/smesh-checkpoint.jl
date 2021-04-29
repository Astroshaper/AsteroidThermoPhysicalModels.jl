

"""
Triangular surface mesh of a polyhedral shape model
"""
struct SMesh
    A::SVector{3,Float64}
    B::SVector{3,Float64}
    C::SVector{3,Float64}
    
    center::SVector{3,Float64}
    normal::SVector{3,Float64}
    area::Float64
    
    f2f::Vector{Int64}
    # ViewFactors::Vector{ViewFactor}
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

    println("Visible faces")
    println(smesh.f2f)
end


SMesh(A, B, C) = SMesh([A, B, C])
SMesh(vs) = SMesh(vs[1], vs[2], vs[3], getcenter(vs), getnormal(vs), getarea(vs), Int64[])

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

getVisibleFaceList(meshes::Vector{SMesh}) = [m.f2f for m in meshes]


################################################################
#                      View Factor
################################################################

struct ViewFactor
    id::Int64     # Index of the interfacing mesh
    fᵢⱼ::Float64  # View factor from mesh i to mesh j
end


"""
    getViewFactor(mᵢ, mⱼ) -> fᵢⱼ

View factor from mesh i to mesh j
assuming Lambertian emission
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


function addViewFactor!(id, meshes, mᵢ, mⱼ)
    
end


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

function isFace(obs::AbstractVector, tar::SMesh)
    (tar.center - obs) ⋅ tar.normal < 0 ? true : false
end


################################################################
#                           Raycast
################################################################

## Viewd from the observing points/mesh
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
    
    PdotE1 = P ⋅ E1
        
    u = (P ⋅ T) / PdotE1
    if 0 ≤ u ≤ 1
        v = (Q ⋅ R) / PdotE1
        if 0 ≤ v ≤ 1
            if 0 ≤ u + v ≤ 1
                t = (Q ⋅ E2) / PdotE1  # 三角形までの距離（t<0なら裏から当たる）
                if t > 0
                    return true
                end
            end
        end
    end
    return false
end


"""
    findVisibleFaces!(obs::SMesh, meshes)

Find faces directly seen from the observer on the face

# Parameters
- `obs`    : Mesh where the observer stands
- `meshes` : Array of SMesh instances
"""
function findVisibleFaces!(obs::SMesh, meshes)
    for i in eachindex(meshes)
        tar = meshes[i]
        isAbove(obs, tar) && isFace(obs, tar) && push!(obs.f2f, i)
    end

    for i in copy(obs.f2f)
        tar_i = meshes[i]
        R = tar_i.center - obs.center
        for j in copy(obs.f2f)
            i == j && continue
            tar_j = meshes[j]

            if raycast(tar_j, R, obs)
                dᵢ = norm(R)                          # distance to i-th mesh
                dⱼ = norm(tar_j.center - obs.center)  # distance to j-th mesh
                dᵢ < dⱼ ? filter!(x->x≠j, obs.f2f) : filter!(x->x≠i, obs.f2f)
            end
        end
    end
end


"""
    findVisibleFaces!(meshes)

Find faces directly seen from the observer on each face

# Parameters
- `meshes` : Array of SMesh instances
"""
function findVisibleFaces!(meshes::Vector{SMesh})
    for obs in meshes
        findVisibleFaces!(obs, meshes)
    end
end


"""
Return true if the observing mesh is illuminated by the direct sunlight, false if not
"""
function isIlluminated(obs::SMesh, r̂☉, meshes::Vector{SMesh})
    for i in obs.f2f
        raycast(meshes[i], r̂☉, obs) && return false
    end
    return true
end


isIlluminated(obs::SMesh, r̂☉, shape) = isIlluminated(obs, r̂☉, shape.smeshes)


function getIlluminatedFaces(r̂☉, meshes::Vector{SMesh})
    illuminated = Int64[]
    for i in eachindex(meshes)
        obs = meshes[i]
        Ψ = obs.normal ⋅ r̂☉  # cosine of the Sun illumination angle
        Ψ > 0 && isIlluminated(obs, r̂☉, meshes) && push!(illuminated, i)
    end
    illuminated
end

function getIlluminatedFaces_PSEUDOCONVEX(r̂☉, meshes::Vector{SMesh})
    illuminated = Int64[]
    for i in eachindex(meshes)
        obs = meshes[i]
        Ψ = obs.normal ⋅ r̂☉  # cosine of the Sun illumination angle
        Ψ > 0 && push!(illuminated, i)
    end
    illuminated
end

