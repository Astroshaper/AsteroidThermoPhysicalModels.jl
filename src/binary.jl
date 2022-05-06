

"""
    mutable struct MutualOrbit{T1}

Mutual orbit parameters of a binary asteroid
"""
mutable struct MutualOrbit{T1, T2}
    a::T1
    e::T1
    I::T1
    # ω::T1
    # Ω::T1
    # Φ::T1
    # tₚ::T1

    # μ::T1
    n::T1
    P::T1

    t::T1
    θ::T1  # Revolutional angle of the mutual orbit
    # u::T1
    # ν::T1
    # r::T2
    # v::T2
    # F☉::T1

    m₁::T1
    m₂::T1
    r₁::T2  # Center of the primary
    r₂::T2  # Center of the secondary
end


"""
    mutable struct Binary{T1, T2, T3, T4}

Describe the state of a binary asteroid

# Fields
- `shape1` : Shape model of the primary
- `shape2` : Shape model of the secondary

- `orbit`        : Orbital elements of the primary
- `mutual_orbit` : Mutual orbit parametes

- `spin1` : Spin parameters of the primary
- `spin2` : Spin parameters of the secondary
"""
mutable struct Binary{T1, T2, T3, T4}
    shape1::T1
    shape2::T1

    orbit::T2
    mutual_orbit::T3

    spin1::T4
    spin2::T4
end


"""
    update_binary!(binary, t)
"""
function update!(binary::Binary, t)
    update_orbit!(binary.orbit, t)

    update_spin!(binary.spin1, t)
    update_spin!(binary.spin2, t)

    binary.mutual_orbit.t = t
    binary.mutual_orbit.θ = binary.mutual_orbit.n * t

    @unpack a, θ, m₁, m₂ = binary.mutual_orbit
    binary.mutual_orbit.r₁ = [cos(θ+π), sin(θ+π), 0] * a * m₂ / (m₁ + m₂)
    binary.mutual_orbit.r₂ = [cos(θ),   sin(θ),   0] * a * m₁ / (m₁ + m₂)
end


"""
    run_binary_TPM!(binary, thermo_params)
"""
function run_binary_TPM!(binary, thermo_params)
    @unpack shape1, shape2, orbit, spin1, spin2 = binary
    @unpack t_bgn, Δt, t_end, P = thermo_params

    init_temps_zero!(shape1, thermo_params)
    init_temps_zero!(shape2, thermo_params)
    
    for t in (t_bgn:Δt:t_end)*P
        update!(binary, t)
        
        r̂☉ = normalize(orbit.r) * -1  # Shift the origin from the sun to the body
        r̂☉₁ = orbit_to_body(r̂☉, spin1)
        r̂☉₂ = orbit_to_body(r̂☉, spin2)
        
        update_flux_sun!(shape1, orbit.F☉, r̂☉₁)
        update_flux_scat_single!(shape1, thermo_params)
        update_flux_rad_single!(shape1, thermo_params)
        
        update_temps!(shape1, thermo_params)
        
        update_flux_sun!(shape2, orbit.F☉, r̂☉₂)
        update_flux_scat_single!(shape2, thermo_params)
        update_flux_rad_single!(shape2, thermo_params)
        
        update_temps!(shape2, thermo_params)
        # println(t, ", ", rad2deg(ϕ₁), ", ", rad2deg(ϕ₂), ", ", r̂☉₁)
    end
end
