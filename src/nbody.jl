

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
    println(io, typeof(p))

    println(" r  : ", p.r)
    println(" v  : ", p.v)
    println(" a  : ", p.a)
    
    println(" a¹ : ", p.a¹)
    println(" a² : ", p.a²)
    println(" a³ : ", p.a³)

    println(" m  : ", p.m)
    println(" R  : ", p.R)
end


function prep_snapshot(ps, params_sim)
    @unpack Δt, t_end, save_interval = params_sim

    num_body = length(ps)
    num_save = length(0:Δt*save_interval:t_end)
    
    ts = Array{Float64}(undef, num_save)    
    rs = Array{Float64}(undef, num_save, 3, num_body)
    vs = Array{Float64}(undef, num_save, 3, num_body)
    as = Array{Float64}(undef, num_save, 3, num_body)
    
    ts, rs, vs, as
end


function save_snapshot!(ts, rs, vs, as, i, t, ps, params_sim)
    @unpack save_interval = params_sim
    (i-1)%save_interval != 0 && return
    
    idx_save = i ÷ save_interval + 1
    ts[idx_save] = t
    for (i, p) in enumerate(ps)
        rs[idx_save, :, i] .= p.r
        vs[idx_save, :, i] .= p.v
        as[idx_save, :, i] .= p.a
    end
end


###################################################################
#                    Barycenter of particles
###################################################################


"""
    getBaryCenter(ps) -> r_G, v_G, a_G

Position, velocity and acceleration of the system's barycenter
"""
function getBaryCenter(ps)
    Σm = sum(ps.m)
    r_G = sum(ps.m .* ps.r) / Σm
    v_G = sum(ps.m .* ps.v) / Σm
    a_G = sum(ps.m .* ps.a) / Σm
    
    r_G, v_G, a_G
end


function setOrigin2BaryCenter!(ps)
    r_G, v_G, a_G = getBaryCenter(ps)
    for p in ps
        p.r .-= r_G
        p.v .-= v_G
        p.a .-= a_G
    end
end


function getBaryCenter(ps, rs, vs, as)
    rs_G = zeros(size(rs[:,:,end]))
    vs_G = zeros(size(vs[:,:,end]))
    as_G = zeros(size(as[:,:,end]))

    for i in 1:size(rs, 3)
        @. rs_G += ps.m[i] * rs[:, :, i]
        @. vs_G += ps.m[i] * vs[:, :, i]
        @. as_G += ps.m[i] * as[:, :, i]
    end

    Σm = sum(ps.m)
    rs_G /= Σm
    vs_G /= Σm
    as_G /= Σm
    
    rs_G, vs_G, as_G
end


function setOrigin2BaryCenter!(ps, rs, vs, as)
    rs_G, vs_G, as_G = getBaryCenter(ps, rs, vs, as)
    
    for i in 1:size(rs, 3)
        @. rs[:,:,i] -= rs_G
        @. vs[:,:,i] -= vs_G
        @. as[:,:,i] -= as_G
    end
end


###################################################################
#                     
###################################################################


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

