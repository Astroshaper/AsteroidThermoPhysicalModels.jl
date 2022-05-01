
"""
    mutable struct Binary{T}

Describe the state of a binary asteroid

# Fields
- `shape1` : Shape model of the primary
- `shape2` : Shape model of the secondary

- `ϕ₁` : Spin phase of the primary
- `ϕ₂` : Spin phase of the secondary
- `ω₁` : Spin velocity of the primary
- `ω₂` : Spin velocity of the secondary
"""
mutable struct Binary{T1, T2, T3}
    shape1::T1
    shape2::T1

    orbit::T2
    # mutual_orbit

    spin1::T3
    spin2::T3
end


function run_binary(shape1, shape2, orbit1, orbit2, spin1, spin2, params_thermo)
    @unpack P, Δt, t_bgn, t_end, Nt, Nz = params_thermo

    init_temperature!(shape1, orbit1, spin1, params_thermo)
    init_temperature!(shape2, orbit2, spin2, params_thermo)

    for t in (t_bgn:Δt:t_end)*P
        ϕ₁ = spin1.ω * t
        ϕ₂ = spin2.ω * t

        F☉, r̂☉ = getSolarCondition(orbit1, spin1, t)
        

        

    end

    # τ̄ = zeros(3)  # Net YORP torque

    # for t in (t_bgn:Δt:t_end)*P
    #     spin_phase = spin.ω * t
    #     F☉, r̂☉ = getSolarCondition(orbit, spin, t)
        
    #     update_flux_sun!(shape, F☉, r̂☉)
    #     update_flux_scat_single!(shape, params_thermo)
    #     update_flux_rad_single!(shape, params_thermo)
        
    #     update_force!(shape, params_thermo)
    #     # update_force_Rubincam!(shape, params_thermo)
        
    #     τ̄ .+= body_to_orbit(SVector{3}(shape.torque), spin.γ, spin.ε, spin_phase)
        
    #     update_temperature!(shape, params_thermo)
    # end
    # τ̄ /= Nt

end