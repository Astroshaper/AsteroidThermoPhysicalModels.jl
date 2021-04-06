## kepler.jl
##     Functions for orbital elements and Kepler motion

using StaticArrays

include("constants.jl")
include("coordinates.jl")


###################################################################
#                       Orbital elements
###################################################################


struct OrbitalElements{T}
    a::T  # semi-major axis
    e::T  # eccentricity
    I::T  # inclination
    ω::T  # argument of periapsis
    Ω::T  # longitude of ascending node
    Φ::T  # mean anomaly   : tₚ = - Φ / n
    tₚ::T  # periapsis passage time

    μ::T  # standard gravitational parameter of the system (GM)
    n::T  # mean motion    : n = √(μ / a^3)
    T::T  # orbital period : T = 2π / n
end


function OrbitalElements(params)

    a = params[:a] * AU      # semi-mojor axis [m]
    e = params[:e]           # eccentricity
    I = deg2rad(params[:I])  # inclination [rad]
    Ω = deg2rad(params[:Ω])  # longitude of the ascending node [rad]
    ω = deg2rad(params[:ω])  # argument of periapsis [rad]
    
    μ = params[:μ]  # Gravitational parameter
    n = √(μ / a^3)   # mean motion  [rad/sec]
    T = 2π / n       # orbital period [sec]
    
    if haskey(params, :Φ)
        Φ  = deg2rad(params[:Φ])  # mean anomaly at the epoch [rad]
        tₚ = - Φ / n              # periapsis passage time [sec]
    elseif haskey(params, :tₚ)
        tₚ = params[:tₚ]
        Φ  = - n * tp
    else
        println("Give [:Φ] or [:tp] in Dict.")
    end

    OrbitalElements(a, e, I, ω, Ω, Φ, tₚ, μ, n, T)
end


function Base.show(io::IO, elms::OrbitalElements)
    println(io, "Orbital elements")
    println(io, "----------------")

    println("Semi-mojor axis        : a  = ", elms.a / AU,     " [AU]")
    println("Eccentricity           : e  = ", elms.e,          " [-]")
    println("Lon. of ascending node : Ω  = ", rad2deg(elms.Ω), " [deg]")
    println("Argument of periapsis  : ω  = ", rad2deg(elms.ω), " [deg]")
    println("Inclination            : I  = ", rad2deg(elms.I), " [deg]")
    println("Periapsis passage time : tₚ = ", elms.tₚ,         " [sec]")
    println("Mean anomaly           : Φ  = ", rad2deg(elms.Φ), " [deg]")

    println()
    println("Other parameters")
    println("----------------")

    println("Gravitational parameter : μ = ", elms.μ,                      " [m^3/s^2]")
    println("Mean motion             : n = ", rad2deg(elms.n) * (3600*24), " [deg/day]")
    println("Orbital period          : T = ", elms.T / (3600*24),          " [day]")
end


###################################################################
#                   Coordinate transformation
#                   R (X, Y, Z) -> r (x, y, z)
###################################################################


ref_to_orb!(r, elms) = ref_to_orb!(r, elms.ω, elms.I, elms.Ω)
ref_to_orb(r, elms) = ref_to_orb(r, elms.ω, elms.I, elms.Ω)


"""
Transform from the reference coordinate system to the orbital plane system

# Parameters
- r : vector to be transformed
- ω : argument of periapsis
- I : inclination
- Ω : longitude of the ascneding node
"""
function ref_to_orb!(r, ω, I, Ω)
    rotateZ!(r, Ω)
    rotateX!(r, I)
    rotateZ!(r, ω)
end


function ref_to_orb(r, ω, I, Ω)
    r = rotateZ(r, Ω)
    r = rotateX(r, I)
    r = rotateZ(r, ω)
end


###################################################################
#                   Coordinate transformation
#                    r(x, y, z) -> R(X, Y, Z)
###################################################################


orb_to_ref!(r, elms) = orb_to_ref!(r, elms.ω, elms.I, elms.Ω)
orb_to_ref(r, elms) = orb_to_ref(r, elms.ω, elms.I, elms.Ω)


"""
Transform from the orbital plane system to the reference coordinate system

# Parameters
- r : vector to be transformed
- ω : argument of periapsis
- I : inclination
- Ω : longitude of the ascneding node
"""
function orb_to_ref!(r, ω, I, Ω)
    rotateZ!(r, -ω)
    rotateX!(r, -I)
    rotateZ!(r, -Ω)
end


function orb_to_ref(r, ω, I, Ω)
    r = rotateZ(r, -ω)
    r = rotateX(r, -I)
    r = rotateZ(r, -Ω)
end


###################################################################
#           R(X, Y, Z), V(Ẋ, Ẏ, Ż) -> Orbital elements
#             @Reference frame
###################################################################


OrbitalElements(R, V, M, m, t) = OrbitalElements(R, V, G*(M+m), t)

function OrbitalElements(R, V, μ, t)

    h = R × V
    
    a = get_a(R, V, μ)
    e = get_e(h, μ, a)
    I = get_I(h)
    Ω = get_Ω(h, I)
    
    ω = get_ω(R, V, a, e, I, Ω)
    f = get_f(R, V, a, e)
    
    E = f2E(f, e)   # eccentric anomaly (E or u)
    n = √(μ / a^3)  # mean motion
    T = 2π / n      # orbital period
    
    tₚ = get_tₚ(t%T, E, e, n)  # time of periapsis passage
    M = - n * tₚ              # mean anomaly (M or Φ)

    OrbitalElements(a, e, I, ω, Ω, M, tₚ, μ, n, T)
end


get_a(R, V, M, m) = get_a(R, V, G*(M+m))
get_a(R, V, μ) = 1 / (2/norm(R) - norm(V)^2/μ)

get_e(R, V, μ, a) = get_e(R × V, μ, a)
get_e(h, μ, a) = √(1 - norm(h)^2/(μ*a))

get_I(R, V) = get_I(R × V)
get_I(h) = acos(h[3]/norm(h))

get_Ω(R, V, I) = get_Ω(R × V, I)

function get_Ω(h, I)
    # sinΩ =  sign(h[3]) * h[1] / (norm(h)*sin(I))
    # cosΩ = -sign(h[3]) * h[2] / (norm(h)*sin(I))
    # Ω = atan(sinΩ, cosΩ)
    
    Ω = atan(sign(h[3])*h[1], -sign(h[3])*h[2])    
    Ω ≥ 0 ? Ω : Ω+2π
end

function get_ω(R, V, a, e, I, Ω)
    ω₊f = get_ω₊f(R, V, I, Ω)
    f = get_f(R, V, a, e)
    
    ω = ω₊f - f
    ω ≥ 0 ? ω : ω+2π
end

function get_ω₊f(R, V, I, Ω)
    I == 0 && return asin((R[2]*cos(Ω) - R[1]*sin(Ω))/norm(R))
    
    sin_ω₊f = R[3] / (norm(R)*sin(I))
    cos_ω₊f = sec(Ω) * (R[1]/norm(R) + sin(Ω)*sin_ω₊f*cos(I))

    ω₊f = atan(sin_ω₊f, cos_ω₊f)
    ω₊f ≥ 0 ? ω₊f : ω₊f + 2π
end

function get_f(R, V, a, e)
    h = R × V
    Ṙ = sign(R ⋅ V) * √(norm(V)^2 - norm(h)^2/norm(R)^2)
    
    sinf = a*(1-e*e) / (norm(h)*e) * Ṙ
    cosf = (a*(1-e*e) / norm(R) - 1) / e
    
    f = atan(sinf, cosf)
    f ≥ 0 ? f : f+2π
end

function get_E(R, a, e)
    E = acos((a - norm(R)) / (a*e))  # eccentric anomaly (E or u)
    E ≥ 0 ? E : E+2π
end

E2f(E, e) = 2 * atan(√((1+e)/(1-e)) * sin(E/2), cos(E/2))
f2E(f, e) = 2 * atan(√((1-e)/(1+e)) * sin(f/2), cos(f/2))

get_tₚ(t, E, e, M, m, a) = get_tₚ(t, E, G*(M+m), a)
get_tₚ(t, E, e, μ, a) = get_tₚ(t, E, e, √(μ/a^3))
get_tₚ(t, E, e, n) = t - (E - e*sin(E)) / n


# ****************************************************************
# ****************************************************************
# ****************************************************************


function getOrbitalElementsError(elms, elms_ref)
    Δa = elms.a - elms_ref.a
    Δe = elms.e - elms_ref.e
    Δω = elms.ω - elms_ref.ω
    ΔT = elms.T - elms_ref.T
    
    return Δa, Δe, Δω, ΔT
end


# ****************************************************************
#           Orbital elements -> r(x, y, z), v(ẋ, ẏ, ż)
#                                @Orbital plane frame
# ****************************************************************


solveKeplerEquation1(elms, t) = solveKeplerEquation1(elms.e, elms.n, elms.tₚ, t)
solveKeplerEquation2(elms, t) = solveKeplerEquation2(elms.e, elms.n, elms.tₚ, t)

get_rv(elms, u) = get_rv(elms.a, elms.e, elms.n, u)
get_r(elms, u) = get_r(elms.a, elms.e, u)
get_v(elms, u) = get_v(elms.a, elms.e, elms.n, u)


"""
    solveKeplerEquation1(e, n, tₚ, t) -> u_i

Solve Kepler's equation with the 1st order method

# Parameters
- `e`  : eccentricity
- `n`  : mean motion
- `tₚ` : periapsis passage time [sec]
- `t`  : time

# Return
- `u_i` : eccentric anomaly (i-th order solution)
"""
function solveKeplerEquation1(e, n, tₚ, t)

    Φ = n * (t - tₚ)  # mean anomaly
    
    u_i = Φ  # eccentirc anomaly at order 0
    
    for i in 1:20
        u_pri = u_i
        u_i = Φ + e * sin(u_i)
        err = abs(1 - u_pri / u_i)
        # println(i, " : ", u_i, " : ", err)
        err < 1e-10 && return u_i
    end
    u_i
end


"""
    solveKeplerEquation2(e, n, tₚ, t) -> u_i

Solve Kepler's equation with the 2nd order method (Newton's method)

# Parameters
- `e`  : eccentricity
- `n`  : mean motion
- `tₚ` : periapsis passage time [sec]
- `t`  : time

# Return
- `u_i` : eccentric anomaly (i-th order solution)
"""
function solveKeplerEquation2(e, n, tₚ, t)

    Φ = n * (t - tₚ)  # mean anomaly
    
    u_i = Φ  # eccentirc anomaly at order 0
    
    for i in 1:20
        u_pri = u_i
        u_i -= (u_i - e * sin(u_i) - Φ) / (1 - e * cos(u_i))
        err = abs(1 - u_pri / u_i)
        # println(i, " : ", u_i, " : ", err)
        err < 1e-10 && return u_i
    end
    u_i
end


"""
    u2ν(u, e) -> ν

Calculate true anomaly ν from eccentric anomaly u and eccentricity e

# Parameters
- `u` : eccentric anomaly
- `e` : eccentricity

# Return
- `ν` : true anomaly

"""
u2ν(u, e) = 2 * atan(√((1+e)/(1-e)) * sin(u*0.5), cos(u*0.5))


"""
    getrv(a, e, n, u) -> r, v

Get a body's positon and velocity based on Kepler's motion

# Parameters
- `a` : semi-major axis
- `e` : eccentricity
- `n` : mean motion
- `u` : eccentric anomaly

# Returns
- `r` : position
- `v` : velocity
(the coordinate system is in the orbital plane and its origin is the eclipse focus)
"""
function get_rv(a, e, n, u)
    sin_u = sin(u)
    cos_u = cos(u)
    
    x = a * (cos_u - e)
    y = a * √(1 - e*e) * sin_u
    z = 0.
    
    r = SA_F64[x, y, z]
    
    x = - a * n * sin_u / (1 - e * cos_u)
    y = a * n * √(1 - e*e) * cos_u / (1 - e * cos_u)
    z = 0.

    v = SA_F64[x, y, z]

    return r, v
end


"""
    get_r(a, e, u) -> r

Get a body's positon based on Kepler's motion

# Parameters
- `a` : semi-major axis
- `e` : eccentricity
- `u` : eccentric anomaly

# Return
- `r` : position
(the coordinate system is in the orbital plane and its origin is the eclipse focus)
"""
function get_r(a, e, u)
    x = a * (cos(u) - e)
    y = a * √(1 - e*e) * sin(u)
    z = 0.

    r = SA_F64[x, y, z]
end

function get_r(orbit, u) = get_r(orbit.a, orbit.e, u)


"""
    get_v(a, e, n, u) -> v

Get a body's velocity based on Kepler's motion

# Parameters
- `a` : semi-major axis
- `e` : eccentricity
- `n` : mean motion
- `u` : eccentric anomaly

# Return
- `v` : velocity
(the coordinate system is in the orbital plane and its origin is the eclipse focus)
"""
function get_v(a, e, n, u)
    sin_u = sin(u)
    cos_u = cos(u)

    x = - a * n * sin_u / (1 - e * cos_u)
    y = a * n * √(1 - e*e) * cos_u / (1 - e * cos_u)
    z = 0.
    
    v = SA_F64[x, y, z]
end

function get_r(orbit, u) = get_r(orbit.a, orbit.e, orbit.n, u)


################################################################
#                     Rotation matrix
################################################################


# Rx(θ) = [
#     1   0      0
#     0  cos(θ) sin(θ)
#     0 -sin(θ) cos(θ)
# ]

# Ry(θ) = [
#     cos(θ) 0 -sin(θ)
#      0     1   0
#     sin(θ) 0  cos(θ)
# ]

# Rz(θ) = [
#      cos(θ) sin(θ) 0
#     -sin(θ) cos(θ) 0
#       0      0     1 
# ]

# orb2ref(r, elms) = orb2ref(r, elms.ω, elms.I, elms.Ω)
# ref2orb(r, elms) = ref2orb(r, elms.ω, elms.I, elms.Ω)

# orb2ref(r, ω, I, Ω) = Rz(-Ω) * Rx(-I) * Rz(-ω) * r
# ref2orb(r, ω, I, Ω) = Rz( ω) * Rx( I) * Rz( Ω) * r


