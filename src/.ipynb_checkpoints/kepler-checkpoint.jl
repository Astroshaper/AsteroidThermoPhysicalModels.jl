## kepler.jl
##     Functions for orbital elements and Kepler motion

using StaticArrays

include("constants.jl")


struct Orbit
    a::Float64   # semi-major axis
    e::Float64   # eccentricity
    I::Float64   # inclination
    ω::Float64   # argument of periapsis
    Ω::Float64   # longitude of ascending node
    Φ::Float64   # mean anomaly   : tp = - Φ / n
    tp::Float64  # periapsis passage time

    μ::Float64  # standard gravitational parameter of the system (GM)
    n::Float64  # mean motion    : n = √(μ / a^3)
    T::Float64  # orbital period : T = 2π / n
end


function Base.show(io::IO, orbit::Orbit)
    # AU = 1.49597870700e11  # [km]

    println(io, "Orbital elements")
    println(io, "----------------")

    println("Semi-mojor axis        : a  = ", orbit.a / AU,     " [AU]")
    println("Eccentricity           : e  = ", orbit.e,          " [-]")
    println("Lon. of ascending node : Ω  = ", rad2deg(orbit.Ω), " [deg]")
    println("Argument of periapsis  : ω  = ", rad2deg(orbit.ω), " [deg]")
    println("Inclination            : I  = ", rad2deg(orbit.I), " [deg]")
    println("Periapsis passage time : tp = ", orbit.tp,         " [sec]")
    println("Mean anomaly           : Φ  = ", rad2deg(orbit.Φ), " [deg]")

    println()
    println("Other parameters")
    println("----------------")

    println("Gravitational parameter : μ = ", orbit.μ,                      " [m^3/s^2]")
    println("Mean motion             : n = ", rad2deg(orbit.n) * (3600*24), " [deg/day]")
    println("Orbital period          : T = ", orbit.T / (3600*24),          " [day]")
end


function setOrbitParams(params::Dict)
    # AU = 1.49597870700e11  # [km]
    # G = 6.6740831e-11      # Gravitational constant [m^3/kg/s^2]
    # GM☉ = 1.32712440018e20

    ## Orbital elements
    a = params[:a] * AU      # semi-mojor axis [km]
    e = params[:e]           # eccentricity
    I = deg2rad(params[:I])  # inclination [rad]
    Ω = deg2rad(params[:Ω])  # longitude of the ascending node [rad]
    ω = deg2rad(params[:ω])  # argument of periapsis [rad]
    
    if haskey(params, :GM)
        GM = params[:GM]
        M  = GM / G
    elseif haskey(params, :M)
        M  = params[:M]
        GM = G * M
    else
        GM = 0.
    end
    
    μ = GM☉ + GM
    n = √(μ / a^3)   # mean motion  [rad/sec]
    T = 2π / n       # orbital period [sec]
    
    if haskey(params, :Φ)
        Φ  = deg2rad(params[:Φ])  # mean anomaly at the epoch [rad]
        tp = - Φ / n              # periapsis passage time [sec]
    elseif haskey(params, :tp)
        tp = params[:tp]
        Φ  = - n * tp
    else
        println("Give [:Φ] or [:tp] in Dict.")
    end

    return Orbit(a, e, I, ω, Ω, Φ, tp, μ, n, T)
end


# ****************************************************************
# ****************************************************************
# ****************************************************************


"""
    solveKeplerEquation1(e, n, tp, t) -> u_i

Solve Kepler's equation with the 1st order method

# Parameters
- `e`  : eccentricity
- `n`  : mean motion
- `tp` : periapsis passage time [sec]
- `t`  : time

# Return
- `u_i` : eccentric anomaly (i-th order solution)
"""
function solveKeplerEquation1(e, n, tp, t)

    Φ = n * (t - tp)  # mean anomaly
    
    u_i = Φ  # eccentirc anomaly at order 0
    
    for i in 1:20
        u_pri = u_i
        u_i = Φ + e * sin(u_i)
        err = abs(1 - u_pri / u_i)
        # println(i, " : ", u_i, " : ", err)

        if err < 1e-10
            return u_i
        end
    end
    u_i
end


"""
    solveKeplerEquation2(e, n, tp, t) -> u_i

Solve Kepler's equation with the 2nd order method (Newton's method)

# Parameters
- `e`  : eccentricity
- `n`  : mean motion
- `tp` : periapsis passage time [sec]
- `t`  : time

# Return
- `u_i` : eccentric anomaly (i-th order solution)
"""
function solveKeplerEquation2(e, n, tp, t)

    Φ = n * (t - tp)  # mean anomaly
    
    u_i = Φ  # eccentirc anomaly at order 0
    
    for i in 1:20
        u_pri = u_i
        u_i -= (u_i - e * sin(u_i) - Φ) / (1 - e * cos(u_i))
        err = abs(1 - u_pri / u_i)
        # println(i, " : ", u_i, " : ", err)

        if err < 1e-10
            return u_i
        end
    end
    u_i
end


solveKeplerEquation1(orbit::Orbit, t) = solveKeplerEquation1(orbit.e, orbit.n, orbit.tp, t)
solveKeplerEquation2(orbit::Orbit, t) = solveKeplerEquation2(orbit.e, orbit.n, orbit.tp, t)


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
function getrv(a, e, n, u)
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
    getposition(a, e, u) -> r

Get a body's positon based on Kepler's motion

# Parameters
- `a` : semi-major axis
- `e` : eccentricity
- `u` : eccentric anomaly

# Return
- `r` : position
(the coordinate system is in the orbital plane and its origin is the eclipse focus)
"""
function getposition(a, e, u)
    x = a * (cos(u) - e)
    y = a * √(1 - e*e) * sin(u)
    z = 0.

    r = SA_F64[x, y, z]
end


"""
    getvelocity(a, e, n, u) -> v

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
function getvelocity(a, e, n, u)
    sin_u = sin(u)
    cos_u = cos(u)

    x = - a * n * sin_u / (1 - e * cos_u)
    y = a * n * √(1 - e*e) * cos_u / (1 - e * cos_u)
    z = 0.
    
    v = SA_F64[x, y, z]
end


getrv(orbit::Orbit, u) = getrv(orbit.a, orbit.e, orbit.n, u)
getposition(orbit::Orbit, u) = getposition(orbit.a, orbit.e, u)
getvelocity(orbit::Orbit, u) = getvelocity(orbit.a, orbit.e, orbit.n, u)




###################################################################
###################################################################
###################################################################


Rx(θ) = [
    1   0      0
    0  cos(θ) sin(θ)
    0 -sin(θ) cos(θ)
]

Ry(θ) = [
    cos(θ) 0 -sin(θ)
     0     1   0
    sin(θ) 0  cos(θ)
]

Rz(θ) = [
     cos(θ) sin(θ) 0
    -sin(θ) cos(θ) 0
      0      0     1 
]

