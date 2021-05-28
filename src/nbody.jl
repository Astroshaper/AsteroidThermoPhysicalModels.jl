

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
    
- `ᵖr`  # Predicted position at next timestep
- `ᵖv`  # Predicted velocity at next timestep

- `ᶜr`  # Corrected position at next timestep
- `ᶜv`  # Corrected position at next timestep

- `m`   # Mass of the particle
- `R`   # Radius of the particle
"""
struct Particle{T1 <: AbstractVector, T2 <: Real}
    r::T1   # Position
    v::T1   # Velocity
    a::T1   # Acceleration
    
    a¹::T1  # 1st dervative of acceleration
    a²::T1  # 2nd dervative of acceleration
    a³::T1  # 3rd dervative of acceleration
    
    ᵖr::T1  # Predicted position at next timestep
    ᵖv::T1  # Predicted velocity at next timestep
    
    ⁺a⁰::T1  # Acceleration at next timestep
    ⁺a¹::T1  # 1st dervative of acceleration at next timestep
    # ⁺a²::T1  # 2nd dervative of acceleration at next timestep
    # ⁺a³::T1  # 3rd dervative of acceleration at next timestep
    
    ᶜr::T1  # Corrected position at next timestep
    ᶜv::T1  # Corrected position at next timestep

    m::T2  # Mass of the particle
    R::T2  # Radius of the particle
end

# struct Particle{T1 <: AbstractVector, T2 <: Real}
#     r::T1  # Position
#     v::T1  # Velocity

#     a⁰::T1  # Acceleration
#     a¹::T1  # 1st dervative of acceleration
#     a²::T1  # 2nd dervative of acceleration
#     a³::T1  # 3rd dervative of acceleration
    
#     ᵖr::T1  # Predicted position at next timestep
#     ᵖv::T1  # Predicted velocity at next timestep
    
#     ⁺a⁰::T1  # Acceleration at next timestep
#     ⁺a¹::T1  # 1st dervative of acceleration at next timestep
#     # ⁺a²::T1  # 2nd dervative of acceleration at next timestep
#     # ⁺a³::T1  # 3rd dervative of acceleration at next timestep
    
#     ᶜr::T1  # Corrected position at next timestep
#     ᶜv::T1  # Corrected position at next timestep

#     m::T2  # Mass of the particle
#     R::T2  # Radius of the particle
# end


function Particle(r, v, m, R)
    
    a   = similar(r) ; a   .= 0
    a¹  = similar(r) ; a¹  .= 0
    a²  = similar(r) ; a²  .= 0
    a³  = similar(r) ; a³  .= 0
    ᵖr  = similar(r) ; ᵖr  .= 0
    ᵖv  = similar(r) ; ᵖv  .= 0
    ⁺a⁰ = similar(r) ; ⁺a⁰ .= 0
    ⁺a¹ = similar(r) ; ⁺a¹ .= 0
    ᶜr  = similar(r) ; ᶜr  .= 0
    ᶜv  = similar(r) ; ᶜv  .= 0
    
    Particle(r, v, a, a¹, a², a³, ᵖr, ᵖv, ⁺a⁰, ⁺a¹, ᶜr, ᶜv, m, R)
end


function Base.show(io::IO, p::Particle)
    println(io, "Particle")

    println(" r  : ", p.r)
    println(" v  : ", p.v)
    println(" a  : ", p.a)
    
    println(" a¹ : ", p.a¹)
    println(" a² : ", p.a²)
    println(" a³ : ", p.a³)
    
    println("ᵖr  : ", p.ᵖr)
    println("ᵖv  : ", p.ᵖv)
    
    println("⁺a⁰ : ", p.⁺a⁰)
    println("⁺a¹ : ", p.⁺a¹)
    
    println("ᶜr  : ", p.ᶜr)
    println("ᶜv  : ", p.ᶜv)

    println("m   : ", p.m)
    println("R   : ", p.R)
end


###################################################################
#                  　　　　 Visualization
###################################################################


get_rs(snapshots, i) = [ps[i].r for ps in snapshots]

get_xs(snapshots, i) = [r[1] for r in get_rs(snapshots, i)]
get_ys(snapshots, i) = [r[2] for r in get_rs(snapshots, i)]
get_zs(snapshots, i) = [r[3] for r in get_rs(snapshots, i)]


"""
Center-of-mass of particles
"""
getParticlesCOM(ps) = sum(ps.m .* ps.r) / sum(ps.m)
