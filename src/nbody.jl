

"""
    Particle{T1 <: AbstractVector, T2 <: Real}

Paticle for N-body simulation

# Fields
- `r`   # Position
- `v`   # Velocity
- `a`   # Acceleration

- `a¹`  # 1st dervative of acceleration
- `a²`  # 2nd dervative of acceleration
- `a³`  # 3rd dervative of acceleration

- `m`   # Mass of the particle
- `R`   # Radius of the particle
"""
struct Particle{T1 <: AbstractVector, T2 <: Real}
    r::T1
    v::T1
    a::T1
    
    a¹::T1
    a²::T1
    a³::T1

    m::T2
    R::T2
end


function Particle(r, v, m, R)
    a   = similar(r) ; a  .= 0.
    a¹  = similar(r) ; a¹ .= 0.
    a²  = similar(r) ; a² .= 0.
    a³  = similar(r) ; a³ .= 0.
    
    Particle(r, v, a, a¹, a², a³, m, R)
end


setParticles(ps...) = StructArray([ps...])


function Base.show(io::IO, p::Particle)
    println(io, "Particle")

    println(" r  : ", p.r)
    println(" v  : ", p.v)
    println(" a  : ", p.a)
    
    println(" a¹ : ", p.a¹)
    println(" a² : ", p.a²)
    println(" a³ : ", p.a³)

    println(" m  : ", p.m)
    println(" R  : ", p.R)
end


function _prep_snapshot(ps, Δt, t_end, save_interval)
    num_body = length(ps)
    num_save = length(0:Δt*save_interval:t_end)
    
    ts = Array{Float64}(undef, num_save)    
    rs = Array{Float64}(undef, num_save, 3, num_body)
    vs = Array{Float64}(undef, num_save, 3, num_body)
    as = Array{Float64}(undef, num_save, 3, num_body)
    
    ts, rs, vs, as
end


function _save_snapshot!(ts, rs, vs, as, i, save_interval, t, ps)
    idx_save = i ÷ save_interval + 1
    ts[idx_save] = t
    for (i, p) in enumerate(ps)
        rs[idx_save, :, i] .= p.r
        vs[idx_save, :, i] .= p.v
        as[idx_save, :, i] .= p.a
    end
end


###################################################################
#                  　　　　 Visualization
###################################################################


get_rs(snapshots, i) = [ps[i].r for ps in snapshots]
get_vs(snapshots, i) = [ps[i].v for ps in snapshots]

get_xs(snapshots, i) = [r[1] for r in get_rs(snapshots, i)]
get_ys(snapshots, i) = [r[2] for r in get_rs(snapshots, i)]
get_zs(snapshots, i) = [r[3] for r in get_rs(snapshots, i)]


"""
"""
function setOrigin!(ps, origin::Particle)
    for p in ps
        p.r .-= origin.r
        p.v .-= origin.v
        p.a .-= origin.a
    end
end


"""
    getBaryCenter(ps) -> r, v, a

Position, velocity and acceleration of the system's barycenter
"""
function getBaryCenter(ps)
    Σm = sum(ps.m)
    r = sum(ps.m .* ps.r) / Σm
    v = sum(ps.m .* ps.v) / Σm
    a = sum(ps.m .* ps.a) / Σm
    
    r, v, a
end


function setOrigin2BaryCenter!(ps)
    r_G, v_G, a_G = getBaryCenter(ps)
    for p in ps
        p.r .-= r_G
        p.v .-= v_G
        p.a .-= a_G
    end
end


getTotalEnergy(ps) = sumKineticEnergy(ps) + sumPotentialEnergy(ps)

getKineticEnergy(p::Particle) = 0.5 * p.m * norm(p.v)^2
sumKineticEnergy(ps) = sum(getKineticEnergy(p) for p in ps)

getPotentialEnergy(pᵢ::Particle, pⱼ::Particle) = - G * pᵢ.m * pⱼ.m / norm(pᵢ.r .- pⱼ.r)

function sumPotentialEnergy(ps)
    E = 0.
    for i in eachindex(ps)
        for j in eachindex(ps)
            i ≥ j && continue
            E += getPotentialEnergy(ps[i], ps[j])
        end
    end
    E
end

