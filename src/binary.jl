

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
# function run_binary_TPM!(binary, thermo_params)
#     @unpack shape1, shape2, orbit, spin1, spin2 = binary
#     @unpack t_bgn, Δt, t_end, P = thermo_params

#     init_temps_zero!(shape1, thermo_params)
#     init_temps_zero!(shape2, thermo_params)

#     ts = (t_bgn:Δt:t_end) * P
#     T₁ = []
#     T₂ = []
    
#     for t in ts
#         update!(binary, t)
        
#         r̂☉ = normalize(orbit.r) * -1  # Shift the origin from the sun to the body
#         r̂☉₁ = orbit_to_body(r̂☉, spin1)
#         r̂☉₂ = orbit_to_body(r̂☉, spin2)
        
#         update_flux_sun!(shape1, orbit.F☉, r̂☉₁)
#         update_flux_scat_single!(shape1, thermo_params)
#         update_flux_rad_single!(shape1, thermo_params)
        
#         update_temps!(shape1, thermo_params)
        
#         update_flux_sun!(shape2, orbit.F☉, r̂☉₂)
#         update_flux_scat_single!(shape2, thermo_params)
#         update_flux_rad_single!(shape2, thermo_params)
        
#         update_temps!(shape2, thermo_params)
#         # println(t, ", ", rad2deg(ϕ₁), ", ", rad2deg(ϕ₂), ", ", r̂☉₁)

#         push!(T₁, surface_temperature(shape1))
#         push!(T₂, surface_temperature(shape2))
#     end
#     ts, T₁, T₂
# end



"""
    run_binary_TPM!(binary, thermo_params)
"""
function run_binary_TPM!(binary, thermo_params)
    @unpack shape1, shape2, orbit, mutual_orbit, spin1, spin2 = binary
    @unpack t_bgn, Δt, t_end, P = thermo_params

    init_temps_zero!(shape1, thermo_params)
    init_temps_zero!(shape2, thermo_params)

    ts = (t_bgn:Δt:t_end) * P
    df = DataFrame(
        t=Float64[], u=Float64[], ν=Float64[],
        f1_x=Float64[], f1_y=Float64[], f1_z=Float64[],
        τ1_x=Float64[], τ1_y=Float64[], τ1_z=Float64[],
        f2_x=Float64[], f2_y=Float64[], f2_z=Float64[],
        τ2_x=Float64[], τ2_y=Float64[], τ2_z=Float64[],
    )
    T₁ = Vector{Float64}[]
    T₂ = Vector{Float64}[]
    
    for t in ts
        update!(binary, t)
        
        r̂☉ = normalize(orbit.r) * -1  # Shift the origin from the sun to the body
        r̂☉₁ = orbit_to_body(r̂☉, spin1)
        r̂☉₂ = orbit_to_body(r̂☉, spin2)
        
        update_flux_sun!(shape1, orbit.F☉, r̂☉₁)
        update_flux_scat_single!(shape1, thermo_params)
        update_flux_rad_single!(shape1, thermo_params)
        
        update_flux_sun!(shape2, orbit.F☉, r̂☉₂)
        update_flux_scat_single!(shape2, thermo_params)
        update_flux_rad_single!(shape2, thermo_params)

        #### Mutual shadowing check ####
        ## Assuming the poles of rotation and mutual orbit coincide

        r̂☉ = rotateZ(r̂☉₁, -spin1.ϕ)  # Sun's direction in the mutable orbit frame
        # r̂☉ = rotateZ(r̂☉₂, -spin2.ϕ)

        r̂₁ = normalize(mutual_orbit.r₁)
        r̂₂ = normalize(mutual_orbit.r₂)
        sep1 = acos(r̂₁ ⋅ r̂☉)
        sep2 = acos(r̂₂ ⋅ r̂☉)
        
        if (t > 100 * spin2.P) && (rad2deg(min(sep1, sep2)) < 20)
            # println(t, ", ", rad2deg(sep1), ", ", rad2deg(sep2))

            for f₁ in shape1.facets
                ## Into the mutable orbit frame
                A₁ = rotateZ(f₁.A, -spin1.ϕ) + mutual_orbit.r₁
                B₁ = rotateZ(f₁.B, -spin1.ϕ) + mutual_orbit.r₁
                C₁ = rotateZ(f₁.C, -spin1.ϕ) + mutual_orbit.r₁
                center1 = rotateZ(f₁.center, -spin1.ϕ) + mutual_orbit.r₁
                # normal1 = rotateZ(f₁.normal, -spin1.ϕ)
                for f₂ in shape2.facets
                    ## Into the mutable orbit frame
                    A₂ = rotateZ(f₂.A, -spin2.ϕ) + mutual_orbit.r₂
                    B₂ = rotateZ(f₂.B, -spin2.ϕ) + mutual_orbit.r₂
                    C₂ = rotateZ(f₂.C, -spin2.ϕ) + mutual_orbit.r₂
                    center2 = rotateZ(f₂.center, -spin2.ϕ) + mutual_orbit.r₂
                    # normal2 = rotateZ(f₂.normal, -spin2.ϕ)

                    raycast(A₂ - center1, B₂ - center1, C₂ - center1, r̂☉) && (f₁.flux.sun = 0)
                    raycast(A₁ - center2, B₁ - center2, C₁ - center2, r̂☉) && (f₂.flux.sun = 0)
                end
            end
        end

        #### Mutual heating check ####
        ## Ignored for now

        update_force!(shape1, thermo_params)
        update_force!(shape2, thermo_params)
        sum_force_torque!(shape1)
        sum_force_torque!(shape2)

        f1 = body_to_orbit(SVector{3}(shape1.force),  spin1)  # Orbital plane frame
        τ1 = body_to_orbit(SVector{3}(shape1.torque), spin1)  # Orbital plane frame

        f2 = body_to_orbit(SVector{3}(shape2.force),  spin2)  # Orbital plane frame
        τ2 = body_to_orbit(SVector{3}(shape2.torque), spin2)  # Orbital plane frame
        
        update_temps!(shape1, thermo_params)
        update_temps!(shape2, thermo_params)

        push!(df, (t, orbit.u, orbit.ν, f1..., τ1..., f2..., τ2...))

        push!(T₁, surface_temperature(shape1))
        push!(T₂, surface_temperature(shape2))
    end
    df, T₁, T₂
end
