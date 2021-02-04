

using LinearAlgebra


struct Particle{T1 <: AbstractVector, T2 <: Real}
    r::T1  # Position
    v::T1  # Velocity

    a⁰::T1  # Acceleration
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


Particle(r, v, m) = Particle(r, v, zeros(3), zeros(3), zeros(3), zeros(3), zeros(3), zeros(3), zeros(3), zeros(3), zeros(3), zeros(3), m, 1.)


Particle(r, v) = Particle(r, v, zeros(3), zeros(3), zeros(3), zeros(3), zeros(3), zeros(3), zeros(3), zeros(3), zeros(3), zeros(3), 1., 1.)


function Base.show(io::IO, p::Particle)
    println(io, "Particle")

    println("r   : ", p.r)
    println("v   : ", p.v)

    println("a⁰  : ", p.a⁰)
    println("a¹  : ", p.a¹)
    println("a²  : ", p.a²)
    println("a³  : ", p.a³)
    
    println("ᵖr  : ", p.ᵖr)
    println("ᵖv  : ", p.ᵖv)
    
    println("⁺a⁰ : ", p.⁺a⁰)
    println("⁺a¹ : ", p.⁺a¹)
    
    println("ᶜr  : ", p.ᶜr)
    println("ᶜv  : ", p.ᶜv)

    println("m   : ", p.m)
    println("R   : ", p.R)
end


#### https://rebound.readthedocs.io/en/latest/c_quickstart.html

# struct reb_particle {
#     double x;           ///< x-position of the particle. 
#     double y;           ///< y-position of the particle. 
#     double z;           ///< z-position of the particle. 
#     double vx;          ///< x-velocity of the particle. 
#     double vy;          ///< y-velocity of the particle. 
#     double vz;          ///< z-velocity of the particle. 
#     double ax;          ///< x-acceleration of the particle. 
#     double ay;          ///< y-acceleration of the particle. 
#     double az;          ///< z-acceleration of the particle. 
#     double m;           ///< Mass of the particle. 
#     double r;           ///< Radius of the particle. 
#     double lastcollision;       ///< Last time the particle had a physical collision.
#     struct reb_treecell* c;     ///< Pointer to the cell the particle is currently in.
#     uint32_t hash;      ///< hash to identify particle.
#     void* ap;           ///< Functionality for externally adding additional properties to particles.
#     struct reb_simulation* sim; ///< Pointer to the parent simulation.
# };

# c.f. StructArrays
# https://github.com/JuliaArrays/StructArrays.jl


###########################################

