

###################################################################
#                         Particles
###################################################################


abstract type AbstractParticle end


"""
    SimpleParticle{T1 <: AbstractVector, T2 <: Real} <: AbstractParticle

Paticle for N-body simulation (e.g. Euler or leapfrog integrator)

# Fields
- `r`   # Position
- `v`   # Velocity
- `a`   # Acceleration

- `m`   # Mass of the particle
- `R`   # Radius of the particle
"""
struct SimpleParticle{T1 <: AbstractVector, T2 <: Real} <: AbstractParticle
    r::T1
    v::T1
    a::T1

    m::T2
    R::T2
end


SimpleParticle(r, v, m, R) = SimpleParticle(r, v, zero(r), m, R)
SimpleParticle(r, v, m) = SimpleParticle(r, v, m, zero(m))
SimpleParticle(p) = SimpleParticle(p.r, p.v, p.a, p.m, p.R)


"""
    Particle{T1 <: AbstractVector, T2 <: Real} <: AbstractParticle

Paticle for N-body simulation (4th-degree Hermite integrator)

# Fields
- `r`   # Position
- `v`   # Velocity
- `a`   # Acceleration

- `a¹`  # 1st dervative of acceleration (jerk)
- `a²`  # 2nd dervative of acceleration
- `a³`  # 3rd dervative of acceleration

- `⁺r`   # Position at next timestep
- `⁺v`   # Velocity at next timestep
- `⁺a`   # Acceleration at next timestep

- `⁺a¹`  # 1st dervative of acceleration at next timestep
- `⁺a²`  # 2nd dervative of acceleration at next timestep
- `⁺a³`  # 3rd dervative of acceleration at next timestep

- `m`   # Mass of the particle
- `R`   # Radius of the particle
"""
struct HermiteParticle{T1 <: AbstractVector, T2 <: Real} <: AbstractParticle
    r::T1
    v::T1
    a::T1
    
    a¹::T1
    a²::T1
    a³::T1
    
    ⁺r::T1
    ⁺v::T1
    ⁺a::T1
    
    ⁺a¹::T1
    ⁺a²::T1
    ⁺a³::T1

    m::T2
    R::T2
end


function HermiteParticle(r, v, a, m, R)
     a¹,  a²,  a³ = zero(r), zero(r), zero(r)
    ⁺r , ⁺v , ⁺a  = zero(r), zero(r), zero(r) 
    ⁺a¹, ⁺a², ⁺a³ = zero(r), zero(r), zero(r)
    
    HermiteParticle(r, v, a, a¹, a², a³, ⁺r, ⁺v, ⁺a, ⁺a¹, ⁺a², ⁺a³, m, R)
end


HermiteParticle(r, v, m, R) = HermiteParticle(r, v, zero(r), m, R)
HermiteParticle(r, v, m) = HermiteParticle(r, v, m, zero(m))
HermiteParticle(p) = HermiteParticle(p.r, p.v, p.a, p.m, p.R)


###################################################################
#                       Initial conditions
###################################################################


setParticles(ps...) = StructArray([ps...])
setParticles(ps::NamedTuple) = StructArray([values(ps)...])


function addParticle!(ps, elms, m, R)
    E = solveKeplerEquation2(elms, rand()*elms.T)
    r, v = get_rv(elms, E)
    r = orb_to_ref(r, elms)
    v = orb_to_ref(v, elms)

    push!(ps, HermiteParticle(collect(r), collect(v), m, R))
end


function initialize_SimpleParticle!(ps, ϵ)
    for p in ps
        p.a .= 0
    end
    
    for i in eachindex(ps)
        for j in eachindex(ps)
            i ≥ j && continue
            ps[j].m == 0 && continue
            
            r = SVector{3}(ps[j].r) - SVector{3}(ps[i].r)
            r² = norm(r)^2 + ϵ^2
            r⁻³ = r²^(-3/2)
            
            ps[i].a .+= G * ps[j].m * r * r⁻³
            ps[j].a .-= G * ps[j].m * r * r⁻³
        end
    end
end


function initialize_HermiteParticle!(ps, ϵ)
    for p in ps
        p.a  .= 0
        p.a¹ .= 0
    end
        
    for i in eachindex(ps)
        for j in eachindex(ps)
            j == i && continue
            ps[j].m == 0 && continue
            
            r = SVector{3}(ps[j].r) - SVector{3}(ps[i].r)
            v = SVector{3}(ps[j].v) - SVector{3}(ps[i].v)
            
            r² = norm(r)^2 + ϵ^2
            r⁻³ = r²^(-3/2)
            r⁻⁵ = r²^(-5/2)
            
            ps[i].a  .+= G * ps[j].m * r * r⁻³
            ps[i].a¹ .+= G * ps[j].m * (v*r⁻³ - 3(v ⋅ r)*r*r⁻⁵)
        end
    end
end


function initialize!(ps, ϵ)
    eltype(ps) <: SimpleParticle && initialize_SimpleParticle!(ps, ϵ)
    eltype(ps) <: HermiteParticle && initialize_HermiteParticle!(ps, ϵ)
end


initialize!(ps) = initialize!(ps, 0)



###################################################################
#                         Simulation
###################################################################


function run_nbody!(ps, params, filepath)
    
    @unpack SOFTENING, TIMESTEP, INTEGRATOR = params
    
    SOFTENING ? (ϵ = params.ϵ) : (ϵ = 0)
    
    initialize!(ps, ϵ)
    
    TIMESTEP == :CONST    && (Δt = params.Δt)
    TIMESTEP == :VARIABLE && (η = params.η)
        
    # INTEGRATOR = :HERMITE4
    # INTEGRATOR = :LEAPFROG
    
    # TIMESTEP = :SHARED
    # TIMESTEP = :INDIVIDUAL
    
    @show SOFTENING
    @show ϵ
    @show TIMESTEP
    @show INTEGRATOR
        
#     open(filepath, "w") do io
        
        
    
#         @unpack η, t_end, save_interval = params_sim

#         ⁺ps = deepcopy(ps)
#         initialize!(ps, ⁺ps)
#         Δt = η * minimum(norm(p.a)/norm(p.a¹) for p in ps)
#         t = 0.
#         i = 0
        
#         E₀ = getTotalEnergy(ps)
#         while t < t_end
#             i%save_interval == 0 && save_snapshot_txt(io, t, ps)
#             collide(ps) != (0, 0) && break
            
#             Hermite4_shared!(ps, ⁺ps, params_sim, Δt)
#             Δt = get_Δt_Aarseth(⁺ps, params_sim)
#             t += Δt
#             i += 1
#         end
#         E₁ = getTotalEnergy(ps)
#         ΔE = (E₁ - E₀)/E₀
#         @show E₀
#         @show E₁
#         @show ΔE
#         println("Collision: ", collide(ps))
#     end
end


###################################################################
#                         Input / Output
###################################################################


function Base.show(io::IO, p::SimpleParticle)
    msg = "$(typeof(p))\n"

    msg *= " r  : $(p.r)\n"
    msg *= " v  : $(p.v)\n"
    msg *= " a  : $(p.a)\n"
    
    msg *= " m  : $(p.m)\n"
    msg *= " R  : $(p.R)\n"
    print(io, msg)
end


function Base.show(io::IO, p::HermiteParticle)
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


"""
    save_snapshot(out::IOStream, t, ps)

# Format of snapshot
    t
    n
    m1 x1 y1 z1 vx1 vy1 vz1
    m2 x2 y2 z2 vx2 vy2 vz2
    :
    :
    mN xN yN zN vxN vyN vzN
"""
function save_snapshot_txt(io::IO, t, ps)
    println(io, t)
    println(io, length(ps))
    for p in ps
        println(io, p.m, ", ",
            p.r[1], ", ", p.r[2], ", ", p.r[3], ", ",
            p.v[1], ", ", p.v[2], ", ", p.v[3]
        )
    end
end


function load_snapshot(filepath)
    ts = Float64[]
    snapshots = []
    open(filepath, "r") do io
        while !eof(io)
            t = parse(Float64, readline(io))
            N = parse(Int64, readline(io))
            ps = StructArray(SimpleParticle{Vector{Float64}, Float64}[])
            for _ in 1:N
                line = parse.(Float64, split(readline(io), ","))
                m = line[1]
                r = line[2:4]
                v = line[5:7]
                p = SimpleParticle(r, v, zeros(3), m, 0.)
                push!(ps, p)
            end
            push!(ts, t)
            push!(snapshots, ps)
        end
    end
    ts, snapshots
end


###################################################################
#                          Collision
###################################################################


distance(p₁, p₂) = norm(SVector{3}(p₁.r) - SVector{3}(p₂.r))

collide(p₁, p₂) = distance(p₁, p₂) < p₁.R + p₂.R ? true : false

function collide(ps)
    for i in eachindex(ps)
        for j in eachindex(ps)
            i ≥ j && continue
            collide(ps[i], ps[j]) && return (i, j)
        end
    end
    return (0, 0)
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

getKineticEnergy(p::AbstractParticle) = 0.5 * p.m * norm(p.v)^2
sumKineticEnergy(ps) = sum(getKineticEnergy(p) for p in ps)

getPotentialEnergy(pᵢ::AbstractParticle, pⱼ::AbstractParticle) = - G * pᵢ.m * pⱼ.m / norm(pᵢ.r .- pⱼ.r)

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

