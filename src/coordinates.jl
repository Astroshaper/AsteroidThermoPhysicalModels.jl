
## coordinate.jl
##     Functions for coordinate transformation


using StaticArrays

# include("kepler.jl")


"""
    eq2ec(α, δ) -> λ, β

Coordinate transform from the equatorial to ecliptic coordinate systems

# Arguments
- `α`: right ascension (RA)
- `δ`: declination (Dec)
# Returns
- `λ`: ecliptic longitude
- `β`: ecliptic latitude
"""
function eq2ec(α, δ)

    ε = deg2rad(23.439279444444445)  # obliquity of the Earth's axis

    sin_ε = sin(ε)
    cos_ε = cos(ε)
        
    sin_α = sin(α)
    cos_α = cos(α)
        
    sin_δ = sin(δ)
    cos_δ = cos(δ)
    
    λ = atan(cos_ε*cos_δ*sin_α + sin_ε*sin_δ, cos_δ*cos_α)
    λ < 0 && (λ += 2π)

    β = asin(-sin_ε*cos_δ*sin_α + cos_ε*sin_δ)

    return λ, β
end


function ec2eq(λ, β)

    ε = deg2rad(23.439279444444445)  # obliquity of the Earth's axis

    sin_ε = sin(ε)
    cos_ε = cos(ε)

    sin_λ = sin(λ)
    cos_λ = cos(λ)
        
    sin_β = sin(β)
    cos_β = cos(β)
    
    α = atan(cos_ε*cos_β*sin_λ - sin_ε*sin_β, cos_β*cos_λ)
    α < 0 && (α += 2π)

    δ = asin(sin_ε*cos_β*sin_λ + cos_ε*sin_β)

    return α, δ
end


# cpdef inline tuple eq2hz(double RA, double Dec, double obslon, double obslat):
#     """
#     Coordinate transform from equatorial to horizontal coordinate system
#         for an observer in (obslon, obslat)
#     """
#     cdef:
#         double hour_angle = obslon - RA  # hour angle = GST + longitude - RA
#         double azimuth, elevation

#     azimuth = atan2(-cos(Dec)*sin(hour_angle), sin(Dec)*cos(obslat) - cos(Dec)*sin(obslat)*cos(hour_angle))

#     if azimuth < 0:
#         azimuth += 2*pi

#     elevation = asin(sin(Dec)*sin(obslat) + cos(Dec)*cos(obslat)*cos(hour_angle))

#     return azimuth, elevation


################################################################
#                 Rotation transformation (inpace)
################################################################

"""
Rotation transformation around X-axis
Parameters:
    v : vector to be transformed
    θ : rotation angle
"""
function rotateX!(v, θ)
    sinθ = sin(θ)
    cosθ = cos(θ)
    
    v[1], v[2], v[3] = v[1], v[2]*cosθ + v[3]*sinθ, - v[2]*sinθ + v[3]*cosθ
end


"""
Rotation transformation around Y-axis
Parameters:
    v : vector to be transformed
    θ : rotation angleon
"""
function rotateY!(v, θ)
    sinθ = sin(θ)
    cosθ = cos(θ)
    
    v[1], v[2], v[3] = v[1]*cosθ - v[3]*sinθ, v[2], v[1]*sinθ + v[3]*cosθ
end


"""
Rotation transformation around Z-axis
Parameters:
    v : vector to be transformed
    θ : rotation angleon
"""
function rotateZ!(v, θ)
    sinθ = sin(θ)
    cosθ = cos(θ)
    
    v[1], v[2], v[3] = v[1]*cosθ + v[2]*sinθ, - v[1]*sinθ + v[2]*cosθ, v[3]
end


################################################################
#                 Rotation transformation (non-inpalce)
################################################################


rotateX(v, θ) = rotateX!(copy(v), θ)
rotateY(v, θ) = rotateY!(copy(v), θ)
rotateZ(v, θ) = rotateZ!(copy(v), θ)


function rotateX(v::SVector, θ)
    sinθ = sin(θ)
    cosθ = cos(θ)
    
    x =   v[1]
    y =   v[2]*cosθ + v[3]*sinθ
    z = - v[2]*sinθ + v[3]*cosθ

    SA_F64[x, y, z]
end


function rotateY(v::SVector, θ)
    sinθ = sin(θ)
    cosθ = cos(θ)
    
    x = v[1]*cosθ - v[3]*sinθ
    y = v[2]
    z = v[1]*sinθ + v[3]*cosθ

    SA_F64[x, y, z]
end


function rotateZ(v::SVector, θ)
    sinθ = sin(θ)
    cosθ = cos(θ)
    
    x =   v[1]*cosθ + v[2]*sinθ
    y = - v[1]*sinθ + v[2]*cosθ
    z =   v[3]

    SA_F64[x, y, z]
end


################################################################
#                   Coordinate transformation
################################################################

"""
Transform the inertial coordinate system to the orbital plane system

# Parameters
- v : vector to be transformed
- ω : argument of periapsis
- I : inclination
- Ω : longitude of the ascneding node
"""
function inertia_to_orbit!(v, ω, I, Ω)
    rotateZ!(v, Ω)
    rotateX!(v, I)
    rotateZ!(v, ω)
end

function inertia_to_orbit(v, ω, I, Ω)
    v = rotateZ(v, Ω)
    v = rotateX(v, I)
    v = rotateZ(v, ω)
end

# inertia_to_orbit!(v, orbit::Orbit) = inertia_to_orbit!(v, orbit.ω, orbit.I, orbit.Ω)
# inertia_to_orbit(v, orbit::Orbit) = inertia_to_orbit(v, orbit.ω, orbit.I, orbit.Ω)


"""
Transform the orbital plane system to the inertial coordinate system

# Parameters
- v : vector to be transformed
- ω : argument of periapsis
- I : inclination
- Ω : longitude of the ascneding node
"""
function orbit_to_inertia!(v, ω, I, Ω)
    rotateZ!(v, -ω)
    rotateX!(v, -I)
    rotateZ!(v, -Ω)
end

function orbit_to_inertia(v, ω, I, Ω)
    v = rotateZ(v, -ω)
    v = rotateX(v, -I)
    v = rotateZ(v, -Ω)
end

# orbit_to_inertia!(v, orbit::Orbit) = orbit_to_inertia!(v, orbit.ω, orbit.I, orbit.Ω)
# orbit_to_inertia(v, orbit::Orbit) = orbit_to_inertia(v, orbit.ω, orbit.I, orbit.Ω)


"""
    orbit_to_body(v::SArray{Tuple{3},Float64,1,3}, γ, ε, ϕ) -> v

# Parameters
- `v` : vector in the orbital plane frame
- `γ` : longitude of vernal equinox direction of the body
- `ε` : obliquity of the spin pole
- `ϕ` : spin phase of the body

# Return
- `v` : vector in the body-fixed frame
"""
function orbit_to_body(v, γ, ε, ϕ)
    v = rotateZ(v, γ)  # body's ecliptic coordinate
    v = rotateX(v, ε)  # body's equatorial coordinate
    v = rotateZ(v, ϕ)  # body-fixed frame
end


function orbit_to_body!(v, γ, ε, ϕ)
    rotateZ!(v, γ)  # body's ecliptic coordinate
    rotateX!(v, ε)  # body's equatorial coordinate
    rotateZ!(v, ϕ)  # body-fixed frame
end


"""
    body_to_orbit(v::SArray{Tuple{3},Float64,1,3}, γ, ε, ϕ) -> v

# Parameters
- `v` : vector in the body-fixed frame
- `γ` : longitude of vernal equinox direction of the body
- `ε` : obliquity of the spin pole
- `ϕ` : spin phase of the body

# Return
- `v` : vector in the orbital plane frame
"""
function body_to_orbit(v, γ, ε, ϕ)
    v = rotateZ(v, -ϕ)  # body's equatorial coordinate
    v = rotateX(v, -ε)  # body's ecliptic coordinate
    v = rotateZ(v, -γ)  # orbital plane frame 
end


function body_to_orbit!(v, γ, ε, ϕ)
    rotateZ!(v, -ϕ)  # body's equatorial coordinate
    rotateX!(v, -ε)  # body's ecliptic coordinate
    rotateZ!(v, -γ)  # orbital plane frame 
end

