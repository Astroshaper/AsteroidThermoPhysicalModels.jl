
# Polyhedral gravitation

using LinearAlgebra
using StaticArrays
using HalfEdges


function setmesh(shapepath)
    topology, _ = loadmesh(shapepath)
    # points *= 1000.
    
    println("Vertices : ", nvertices(topology))
    println("Faces    : ", nfaces(topology))
    println("Edges    : ", nedges(topology))
end


"""
    getFaceCenters(topology, points) -> facecenters

Get an array of center-positions of all faces

# Arguments
- `topology` : Topology instace of HalfEdges.jl
- `points`   : vertex-list of the polyhedron
# Return
- facecenters : array of center-positions of faces
"""
function getFaceCenters(topology, points)
    num_face = nfaces(topology)
    faces = facelist(topology)
    facecenters = Array{SVector{3,Float64}, 1}(undef, num_face)

    for i in 1:num_face
        face = faces[i]
        center = points[face[1]] + points[face[2]] + points[face[3]]
        center /= 3
        facecenters[i] = center
    end
    return facecenters
end


"""
    getFaceAreas(topology, points) -> faceareas

Get an array of areas of all faces

# Arguments
- `topology` : Topology instace of HalfEdges.jl
- `points`   : vertex-list of the polyhedron
# Return
- faceareas : array of areas of faces
"""
function getFaceAreas(topology, points)
    num_face = nfaces(topology)
    faceareas = Array{Float64, 1}(undef, num_face)
    
    for i in 1:num_face
        area = HalfEdges.area(topology, points, FaceHandle(i))
        faceareas[i] = area
    end
    return faceareas
end


"""
    getFaceNormals(topology, points) -> facenormals

Get an array of normal vectors to all faces

# Arguments
- `topology` : Topology instace of HalfEdges.jl
- `points`   : vertex-list of the polyhedron
# Return
- facenormals : array of face normals
"""
getFaceNormals(topology, points) = HalfEdges.normals(topology, points)


"""
    getFaceVertices(topology, points) -> r1s, r2s, r3s
"""
function getFaceVertices(topology, points)
    num_face = nfaces(topology)
    faces = facelist(topology)
    
    r1s = Array{SVector{3,Float64}, 1}(undef, num_face)
    r2s = Array{SVector{3,Float64}, 1}(undef, num_face)
    r3s = Array{SVector{3,Float64}, 1}(undef, num_face)

    for i in 1:num_face
        r1s[i] = points[faces[i][1]]
        r2s[i] = points[faces[i][2]]
        r3s[i] = points[faces[i][3]]
    end
    return r1s, r2s, r3s
end


"""
    getHalfEdgeVertices(topology, points) -> heads, tails
"""
function getHalfEdgeVertices(topology, points)
    num_halfedge = nhalfedges(topology)
    
    heads = Array{SVector{3,Float64}, 1}(undef, num_halfedge)
    tails = Array{SVector{3,Float64}, 1}(undef, num_halfedge)
    
    for i in 1:num_halfedge
        he = topology.he[i]

        heads[i] = points[he.head]                       # head of the halfedge
        tails[i] = points[HalfEdges.tail(topology, he)]  # tail of the halfedge
    end
    return heads, tails
end



################################################################
#                Constant density polyhedron
################################################################


getSurfaceGravity(topology, points, ρ) = getgravity(topology, points, ρ, getFaceCenters(topology, points))


function getgravity(topology, points, ρ, positions)
    num_halfedge = nhalfedges(topology)
    num_face = nfaces(topology)
    num_data = length(positions)

    potentials = zeros(num_data)
    accelerations = [SA_F64[0, 0, 0] for i in 1:num_data]
    # accelerations = [[0., 0., 0.] for i in 1:num_data]

    heads, tails = getHalfEdgeVertices(topology, points)
    halfedgedyads = getHalfEdgeDyads(topology, points)
    
    r1s, r2s, r3s = getFaceVertices(topology, points)
    facedyads = getFaceDyads(topology, points)
    #### ここまで 118.81 k allocations: 6.495 MiB ####
    
    for i in eachindex(positions)
        pos = positions[i]
        
        for j in 1:num_halfedge
            r_i = heads[j]
            r_j = tails[j]
            E_he = halfedgedyads[j]

            pot, acc = getHalfEdgeTerm(r_i, r_j, E_he, pos)
            potentials[i]    += pot
            accelerations[i] += acc
        end

        for j in 1:num_face
            r_i = r1s[j]
            r_j = r2s[j]
            r_k = r3s[j]
            F_f = facedyads[j]

            pot, acc = getFaceTerm(r_i, r_j, r_k, F_f, pos)
            potentials[i]    -= pot
            accelerations[i] -= acc
        end
    end

    G = 6.6740831e-11  # Gravitational constant [m^3/kg/s^2]
    potentials    *= -0.5 * G * ρ
    accelerations *=        G * ρ

    return potentials, accelerations
end


function getHalfEdgeTerm(r_i, r_j, E_he, pos)
    r_he = r_i - pos
    L_he = getEdgeFactor(r_i, r_j, pos)

    acc = E_he * r_he .* L_he  # dyad, vectorの積でallocationが発生する
    pot = r_he ⋅ acc

    return pot, acc
end


function getFaceTerm(r_i, r_j, r_k, F_f, pos)
    r_f = r_i - pos
    ω_f = getFaceFactor(r_i, r_j, r_k, pos)
    
    acc = F_f * r_f .* ω_f  # vector, dyadの積でallocationが発生する
    pot = r_f ⋅ acc

    return pot, acc
end


################################################################
################################################################
################################################################


"""
"""
# function getSurfaceGravity(topology, points, ρ)
#     num_face = nfaces(topology)

#     facecenters = getFaceCenters(topology, points)
#     halfedgedyads = getHalfEdgeDyads(topology, points)
#     facedyads = getFaceDyads(topology, points)

#     potential = zeros(num_face)
#     attraction = zeros(num_face, 3)

#     for i in 1:num_face
#         U, g = getGravity(topology, points, halfedgedyads, facedyads, facecenters[i], ρ)
#         potential[i] = U
#         attraction[i, :] = g
#     end

#     return potential, attraction
# end


"""
    getGravity(topology, points, halfedgedyads, facedyads, pos, ρ=1.) -> U, g

Calculate the gravity field of a constant density polyhedron

# Arguments
- `topology`      : Topology instace of HalfEdges.jl
- `points`        : vertex-list of the polyhedron
- `halfedgedyads` : array of halfedge dyads, E_he
- `facedyads`     : array of face dyads, F_f
- `pos`           : position where you want the gravity
- `ρ`             : constant density of the polyhedron (default: 1.0)
# Returns
- `U` : potential
- `g` : attraction
"""
# function getGravity(topology, points, halfedgedyads, facedyads, pos, ρ)
    
#     U_e, g_e = sumHalfEdgeTerms(getEdgeVertices(topology, points)..., halfedgedyads, pos)

#     U_f, g_f = sumFaceTerms(getFaceVertices(topology, points)..., facedyads, pos)

#     G = 6.6740831e-11  # Gravitational constant [m^3/kg/s^2]
#     U = -0.5 * G * ρ * (U_e - U_f)
#     g =        G * ρ * (g_e - g_f)

#     return U, g
# end


"""
    sumHalfEdgeTerms(topology, points, halfedgedyads, pos) -> U_he, g_he

Get edge terms of potential and attraction

# Aruguments
- `topology`  : Topology instace of HalfEdges.jl
- `points`        : vertex-list of the polyhedron
- `halfedgedyads` : array of halfedge dyads, E_he
- `pos`           : position where you want the gravity
# Returns
- `U_he` : halfedge term of potential
- `g_he` : halfedge term of attraction
"""
# function sumHalfEdgeTerms(heads, tails, halfedgedyads, pos)
#     num_halfedge = size(halfedgedyads, 1)

#     U_he = 0.
#     g_he = [0., 0., 0.]

#     for i in 1:num_halfedge
#         r_i = heads[i]
#         r_j = tails[i]
#         E_he = halfedgedyads[i]

#         r_he = r_i .- pos
#         L_he = getEdgeFactor(r_i, r_j, pos)

#         U_he += r_he ⋅ (E_he * r_he) * L_he
#         g_he +=        (E_he * r_he) * L_he
#     end
#     return U_he, g_he
# end


"""
    sumFaceTerms(topology, points, facedyads, pos) -> U_f, g_f

Sum face terms over all faces

# Aruguments
- `points`    : vertex-list of the polyhedron
- `faces`     : list of face vertices
- `facedyads` : array of face dyads, F_f
- `pos`       : position where you want the gravity
# Returns
- `U_f` : face term of potential
- `g_f` : face term of attraction
"""
# function sumFaceTerms(r1s, r2s, r3s, facedyads, pos)
#     num_face = size(facedyads, 1)

#     U_f = 0.
#     g_f = [0., 0., 0.]
    
#     for i in 1:num_face
#         r_i = r1s[i]
#         r_j = r2s[i]
#         r_k = r3s[i]

#         F_f = facedyads[i]

#         r_f = r_i .- pos
#         ω_f = getFaceFactor(r_i, r_j, r_k, pos)

#         U_f += r_f ⋅ (F_f * r_f) * ω_f
#         g_f +=       (F_f * r_f) * ω_f
#     end
#     return U_f, g_f
# end


"""
    getEdgeFactor(r_i, r_j) -> L_e

Calculate a edge factor, L_e

# Arguments
- `r_i` : start-point of the edge
- `r_j` : end-point of the edge
- `pos` : position where you want the gravity (if skipped, the origin is given)
# Returns
- `L_e`: edge factor
"""
function getEdgeFactor(r_i, r_j)
    r_i_norm = norm(r_i)
    r_j_norm = norm(r_j)
    e_ij_norm = norm(r_i .- r_j)

    L_e = r_i_norm + r_j_norm + e_ij_norm
    L_e /= r_i_norm + r_j_norm - e_ij_norm
    L_e = log(L_e)
end

getEdgeFactor(r_i, r_j, pos) = getEdgeFactor(r_i .- pos, r_j .- pos)


"""
    getFaceFactor(r_i, r_j, r_k) -> ω_f

Calculate a face factor, ω_f

# Arguments
- `r_i` : 1st vertex of the face (counter-clockwise)
- `r_j` : 2nd vertex of the face
- `r_k` : 3rd vertex of the face 
- `pos` : position where you want the gravity (if skipped, the origin is given)
# Returns
- `ω_f`: face factor
"""
function getFaceFactor(r_i, r_j, r_k)
    r_i_norm = norm(r_i)
    r_j_norm = norm(r_j)
    r_k_norm = norm(r_k)

    y = r_i ⋅ (r_j × r_k)

    x = r_i_norm * r_j_norm * r_k_norm
    x += r_i_norm * (r_j ⋅ r_k)
    x += r_j_norm * (r_k ⋅ r_i)
    x += r_k_norm * (r_i ⋅ r_j)

    ω_f = 2 * atan(y, x)
end

getFaceFactor(r_i, r_j, r_k, pos) = getFaceFactor(r_i .- pos, r_j .- pos, r_k .- pos)


"""
    getHalfEdgeDyads(topology, points) -> halfedgedyads

Get an array of dyads for all halfedges

# Arguments
- `topology` : Topology instace of HalfEdges.jl
- `points`   : vertex-list of the polyhedron
# Return
- halfedgedyads : array of halfedge dyads
"""
function getHalfEdgeDyads(topology, points)
    num_halfedge = nhalfedges(topology)
    halfedgedyads = Array{SMatrix{3,3,Float64}, 1}(undef, num_halfedge)

    for i in 1:num_halfedge
        halfedge = topology.he[i]

        _head = points[halfedge.head]
        _tail = points[HalfEdges.tail(topology, halfedge)]
        edge = normalize(_head - _tail)  # normalized vector from tail to head

        facenormal = HalfEdges.trinormal(topology, points, HalfEdgeHandle(i))
        edgenormal = edge × facenormal
        E_he = getdyad(facenormal, edgenormal)
        halfedgedyads[i] = E_he
    end
    return halfedgedyads
end


"""
    getFaceDyads(topology, points) -> facedyads

Get an array of dyads for all faces

# Arguments
- `topology` : Topology instace of HalfEdges.jl
- `points`   : vertex-list of the polyhedron
# Return
- facedyads : array of face dyads
"""
function getFaceDyads(topology, points)
    num_face = nfaces(topology)
    facedyads = Array{SMatrix{3,3,Float64}, 1}(undef, num_face)
    facenormals = HalfEdges.normals(topology, points)

    for i in 1:num_face
        n̂ = facenormals[i]
        F_f = getdyad(n̂, n̂)
        facedyads[i] = F_f
    end
    return facedyads
end


"""
    getdyad(a, b) -> dyad

Make a dyad from two vectors.
You can simply write: `a * b'`

# Arguments
- `a`, `b`: two vectors
# Return
- `dyad`: dyad
"""
function getdyad(a, b)
    dyad = SA_F64[
        a[1]*b[1] a[1]*b[2] a[1]*b[3];
        a[2]*b[1] a[2]*b[2] a[2]*b[3];
        a[3]*b[1] a[3]*b[2] a[3]*b[3];
    ]
end