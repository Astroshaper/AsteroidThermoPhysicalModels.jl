
###################################################################
#                       Orbital elements
###################################################################

"""
# Orbital elements
- `a`  : Semi-major axis
- `e`  : Eccentricity
- `I`  : Inclination
- `ω`  : Argument of periapsis
- `Ω`  : Longitude of ascending node
- `Φ`  : Mean anomaly (tₚ = - Φ / n)
- `tₚ`  : Periapsis passage time

- `μ`  : Standard gravitational parameter of the system (GM)
- `n`  : Mean motion    : n = √(μ / a^3)
- `P`  : Orbital period : T = 2π / n

## Time-variables
- `t`  : Time
- `u`  : Eccentric anomaly
- `ν`  : True anomaly
- `r`  : Position in the orbital plane frame (located at the focus of the ellipse)
- `v`  : Velocity in the orbital plane frame (located at the focus of the ellipse)
- `F☉` : Solar irradiation at the orbit [W/m²]
"""
mutable struct OrbitalElements{T1, T2}
    a::T1
    e::T1
    I::T1
    ω::T1
    Ω::T1
    Φ::T1
    tₚ::T1

    μ::T1
    n::T1
    P::T1

    t::T1
    u::T1
    ν::T1
    r::T2
    v::T2
    F☉::T1
end


function OrbitalElements(params; t=0.)

    a = params[:a] * AU      # Semi-mojor axis [m]
    e = params[:e]           # Eccentricity
    I = deg2rad(params[:I])  # Inclination [rad]
    Ω = deg2rad(params[:Ω])  # Longitude of the ascending node [rad]
    ω = deg2rad(params[:ω])  # Argument of periapsis [rad]
    
    μ = params[:μ]  # Gravitational parameter
    n = √(μ / a^3)  # Mean motion [rad/sec]
    P = 2π / n      # Orbital period [sec]
    
    if haskey(params, :Φ)
        Φ  = deg2rad(params[:Φ])  # Mean anomaly at the epoch [rad]
        tₚ = - Φ / n              # Periapsis passage time [sec]
    elseif haskey(params, :tₚ)
        tₚ = params[:tₚ]
        Φ  = - n * tₚ
    else
        println("Give [:Φ] or [:tp] in Dict.")
    end

    u = solveKeplerEquation2(e, n, tₚ, t)
    ν = u2ν(u, e)
    r, v = get_rv(a, e, n, u)
    F☉ = SOLAR_CONST / (norm(r) / AU)^2

    OrbitalElements(a, e, I, ω, Ω, Φ, tₚ, μ, n, P, t, u, ν, r, v, F☉)
end


function OrbitalElements(a, e, I, μ)
    a *= AU
    I = deg2rad(I)
    ω = 2π * rand()
    Ω = 2π * rand()
    Φ = 2π * rand()
    
    n = √(μ / a^3)
    P = 2π / n
    tₚ = - Φ / n

    t = rand(0:T)
    u = solveKeplerEquation2(e, n, tₚ, t)
    ν = u2ν(u, e)
    r, v = get_rv(a, e, n, u)
    F☉ = SOLAR_CONST / (norm(r) / AU)^2

    OrbitalElements(a, e, I, ω, Ω, Φ, tₚ, μ, n, P, t, u, ν, r, v, F☉)
end


function Base.show(io::IO, orbit::OrbitalElements)
    @unpack a, e, I, ω, Ω, Φ, tₚ, μ, n, P, t, u, ν, r, v, F☉ = orbit

    msg = "--------------------\n"
    msg *= "  Orbital elements  \n"
    msg *= "--------------------\n"

    msg *= "    Semi-mojor axis         : a  = $(a / AU    ) [AU]\n"
    msg *= "    Eccentricity            : e  = $(e         ) [-]\n"
    msg *= "    Lon. of ascending node  : Ω  = $(rad2deg(Ω)) [deg]\n"
    msg *= "    Argument of periapsis   : ω  = $(rad2deg(ω)) [deg]\n"
    msg *= "    Inclination             : I  = $(rad2deg(I)) [deg]\n"
    msg *= "    Periapsis passage time  : tₚ = $(tₚ        ) [sec]\n"
    msg *= "    Mean anomaly            : Φ  = $(rad2deg(Φ)) [deg]\n"

    msg *= "--------------------\n"
    msg *= "  Other parameters  \n"
    msg *= "--------------------\n"

    msg *= "    Gravitational parameter : μ = $(μ                     ) [m^3/s^2]\n"
    msg *= "    Mean motion             : n = $(rad2deg(n) * (3600*24)) [deg/day]\n"
    msg *= "    Orbital period          : P = $(P / (3600*24)         ) [day]\n"

    msg *= "------------------\n"
    msg *= "  Time-variables  \n"
    msg *= "------------------\n"
    msg *= "    Time                    : t  = $(t         ) [sec]\n"
    msg *= "    Eccentric anomaly       : u  = $(rad2deg(u)) [deg]\n"
    msg *= "    True anomaly            : ν  = $(rad2deg(ν)) [deg]\n"
    msg *= "    Position                : r  = $(r         ) [m]\n"
    msg *= "    Velocity                : v  = $(v         ) [m/s]\n"
    msg *= "    Solar irradiation       : F☉ = $(F☉       ) [W/m²]\n"
    print(io, msg)
end


###################################################################
#                   Update orbital parameters
###################################################################

"""
    update_orbit!(orbit::OrbitalElements, t)
"""
function update_orbit!(orbit::OrbitalElements, t)
    u = solveKeplerEquation2(orbit, t)
    ν = u2ν(u, orbit)
    r, v = get_rv(orbit, u)
    F☉ = SOLAR_CONST / (norm(r) / AU)^2

    orbit.t  = t
    orbit.u  = u
    orbit.ν  = ν
    orbit.r  = r
    orbit.v  = v
    orbit.F☉ = F☉
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


function OrbitalElements(R::AbstractVector, V::AbstractVector, μ, t)

    h = R × V
    
    a = get_a(R, V, μ)
    e = get_e(h, μ, a)
    I = get_I(h)
    Ω = get_Ω(h, I)
    
    ω = get_ω(R, V, a, e, I, Ω)
    f = get_f(R, V, a, e)
    
    E = f2E(f, e)   # eccentric anomaly (E or u)
    n = √(μ / a^3)  # mean motion
    P = 2π / n      # orbital period
    
    tₚ = get_tₚ(t%T, E, e, n)  # time of periapsis passage
    M = - n * tₚ              # mean anomaly (M or Φ)

    OrbitalElements(a, e, I, ω, Ω, M, tₚ, μ, n, P)
end


function OrbitalElements(p1, p2, t)
    ps = setParticles(p1, p2)
    r_G, v_G, a_G = getBaryCenter(ps)
    OrbitalElements(p2.r - r_G, p2.v - v_G, p1.m, p2.m, t)
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
    # ω ≥ 0 ? ω : ω+2π
end

function get_ω₊f(R, V, I, Ω)
    I == 0 && return asin((R[2]*cos(Ω) - R[1]*sin(Ω))/norm(R))
    
    sin_ω₊f = R[3] / (norm(R)*sin(I))
    cos_ω₊f = sec(Ω) * (R[1]/norm(R) + sin(Ω)*sin_ω₊f*cos(I))

    ω₊f = atan(sin_ω₊f, cos_ω₊f)
    # ω₊f ≥ 0 ? ω₊f : ω₊f + 2π
end

function get_f(R, V, a, e)
    h = R × V
    Ṙ = sign(R ⋅ V) * √(norm(V)^2 - norm(h)^2/norm(R)^2)
    
    sinf = a*(1-e*e) / (norm(h)*e) * Ṙ
    cosf = (a*(1-e*e) / norm(R) - 1) / e
    
    f = atan(sinf, cosf)
    # f ≥ 0 ? f : f+2π
end

function get_E(R, a, e)
    E = acos((a - norm(R)) / (a*e))  # eccentric anomaly (E or u)
    # E ≥ 0 ? E : E+2π
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


solveKeplerEquation1(elms, ts::T) where T<:AbstractArray = [solveKeplerEquation1(elms, t) for t in ts]
solveKeplerEquation2(elms, ts::T) where T<:AbstractArray = [solveKeplerEquation2(elms, t) for t in ts]

solveKeplerEquation1(elms, t) = solveKeplerEquation1(elms.e, elms.n, elms.tₚ, t)
solveKeplerEquation2(elms, t) = solveKeplerEquation2(elms.e, elms.n, elms.tₚ, t)

get_rv(elms, u) = get_rv(elms.a, elms.e, elms.n, u)
get_r(elms, u) = get_r(elms.a, elms.e, u)
get_v(elms, u) = get_v(elms.a, elms.e, elms.n, u)

heliocentric_distance(elms, u) = norm(get_r(elms, u))
heliocentric_distance(elms, us::T) where T<:AbstractArray = [heliocentric_distance(elms, u) for u in us]


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
- `u` : Eccentric anomaly
- `e` : Eccentricity

# Return
- `ν` : True anomaly

"""
u2ν(u, e) = 2 * atan(√((1+e)/(1-e)) * sin(u*0.5), cos(u*0.5))
u2ν(u, elms::OrbitalElements) = u2ν(u, elms.e)

u2ν(us::T, e) where T<:AbstractArray = [u2ν(u, e) for u in us]
u2ν(us::T, elms::OrbitalElements) where T<:AbstractArray = u2ν(us, elms.e)


"""
    get_rv(a, e, n, u) -> r, v

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


