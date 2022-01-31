

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