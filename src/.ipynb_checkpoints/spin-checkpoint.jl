

"""
    Spin{T1, T2}

Spin parameters of an asteroid

# Fields

## Spin pole @ equatorial coordinate
- `α`  # Right ascension (RA)
- `δ`  # Declination (Dec)

## Spin pole @ ecliptic coordinate
- `λ`  # Ecliptic longitude
- `β`  # Ecliptic latitude

## Other parameters
- `T`  # Spin period [sec]
- `ω`  # Angular velocity [rad/sec]
- `ŝ`  # Spin pole direction (normalized)
- `ε`  # Obliquity
- `γ`  # Vernal equinox lon. from the direction of perihelion
"""
struct Spin{T1, T2}
    α::T1
    δ::T1

    λ::T1
    β::T1

    T::T1
    ω::T1
    ŝ::T2
    ε::T1

    γ::T1
end


function Base.show(io::IO, spin::Spin)
    println(io, "Spin parameters")
    println(io, "---------------")

    println("Right ascension (RA) : α = ", rad2deg(spin.α), " [deg]")
    println("Declination (Dec)    : δ = ", rad2deg(spin.δ), " [deg]")
    println("Ecliptic longitude   : λ = ", rad2deg(spin.λ), " [deg]")
    println("Ecliptic latitude    : β = ", rad2deg(spin.β), " [deg]")
    println("Obliquity            : ε = ", rad2deg(spin.ε), " [deg]")
    println("Spin period          : P = ", spin.T / 3600,   " [hours]")
    println("Spin rate            : ω = ", spin.ω,          " [rad/sec]")
    println("Vernal equinox lon.  : γ = ", rad2deg(spin.γ), " [deg]")
    println("                           (longitude from the periheion direction)")
end


function setSpinParams(params, orbit)

    if haskey(params, :α) && haskey(params, :δ)
        α = deg2rad(params[:α])
        δ = deg2rad(params[:δ])
        λ, β = eq2ec(α, δ)
    elseif haskey(params, :λ) && haskey(params, :β)
        λ = deg2rad(params[:λ])
        β = deg2rad(params[:β])
        α, δ = ec2eq(λ, β)
    else
        println("Give spin pole direction.")
    end
    
    ŝ = getSpinNormal(λ, β, orbit)  # Orbital plane frame
    ε = acos(ŝ[3])
    γ = getVernalEquinox(ŝ)

    T = params[:T] * 3600
    ω = 2π / T

    return Spin(α, δ, λ, β, T, ω, ŝ, ε, γ)
end


function setSpinParams(T::Float64, ŝ::AbstractVector, orbit)
    
    T *= 3600
    ω = 2π / T
    ε = acos(ŝ[3])
    γ = getVernalEquinox(ŝ)
    
    λ = - atan(ŝ[2], ŝ[1])
    β = asin(ŝ[3])
    α, δ = ec2eq(λ, β)

    return Spin(α, δ, λ, β, T, ω, ŝ, ε, γ) 
end


##############################################################################
##############################################################################


"""
Get a spin pole direction in the orbital plane frame

# Paramters
- λ : ecliptic longitude
- β : ecliptic latitude
- orbit : orbital elements
# Return
- ŝ : spin pole direction (normalized)
"""
function getSpinNormal(λ, β, orbit)
    ŝ = SA_F64[cos(β) * cos(λ), - cos(β) * sin(λ), sin(β)]  # inertial frame
    ŝ = inertia_to_orbit(ŝ, orbit)                          # orbital plane frame
end


"""
Get a unit vector ê⟂1 and ê⟂2
"""
function getSpinUnits(spin::Spin)
    N̂ = SA_F64[0, 0, 1]
    ê1 = (spin.ŝ * cos(spin.ε) - N̂) / sin(spin.ε)
    ê2 = (spin.ŝ × N̂) / sin(spin.ε)
    
    ê1, ê2
end


"""
Get a unit vector ê⟂1 only
"""
function getSpinUnit1(spin::Spin)
    N̂ = SA_F64[0, 0, 1]
    ê1 = (spin.ŝ * cos(spin.ε) - N̂) / sin(spin.ε)
end


"""
Get a unit vector ê⟂2 only
"""
function getSpinUnit2(spin::Spin)
    N̂ = SA_F64[0, 0, 1]
    ê2 = (spin.ŝ × N̂) / sin(spin.ε)
end


"""
Get a longitude of the vernal equinox with respect to the perihelion
"""
function getVernalEquinox(spin::Spin)    
    ê2 = getSpinUnit2(spin)
    vernal_equinox = atan(ê2[2], ê2[1]) + π
end

function getVernalEquinox(ŝ)
    N̂ = SA_F64[0, 0, 1]
    ê2 = ŝ × N̂
    vernal_equinox = atan(ê2[2], ê2[1]) + π
end


"""
Get a longitude of the autumnal equinox with respect to the perihelion
"""
function getAutumnalEquinox(spin)
    ê2 = getSpinUnit2(spin)
    autumnal_equinox = atan(ê2[2], ê2[1])
    autumnal_equinox < 0 && (autumnal_equinox += 2π)
    autumnal_equinox
end
