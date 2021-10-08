

"""
    run_leapfrog(ps, params_sim, filename)
"""
function run_leapfrog(ps, params_sim, filename)
    open(filename, "w") do out
    
        @unpack ϵ, Δt, t_end, Δt_out = params_sim
        calc_force!(ps, ϵ)

        t = 0.
        E₀ = getTotalEnergy(ps)
        while t < t_end
            t%Δt_out == 0 && save_snapshot_txt(out, t, ps)
            
            leapfrog!(ps, ϵ, Δt)
            t += Δt
        end
        E₁ = getTotalEnergy(ps)
        ΔE = E₁ - E₀
        @show E₀
        @show E₁
        @show ΔE / E₀
    end
end


function calc_force!(ps, ϵ)
    for p in ps
        p.a .= 0
    end
    
    for i in eachindex(ps)
        for j in eachindex(ps)
            i ≥ j && continue
            
            r = SVector{3}(ps[j].r) - SVector{3}(ps[i].r)
            r² = norm(r)^2 + ϵ^2
            r⁻³ = r²^(-3/2)
            
            ps[i].a .+= G * ps[j].m * r * r⁻³
            ps[j].a .-= G * ps[i].m * r * r⁻³
        end
    end
end


function leapfrog!(ps, ϵ, Δt)
    for p in ps
        @. p.v += p.a * 0.5Δt  # Velocity at half-step later
        @. p.r += p.v * Δt     # Position at one-step later
    end
    
    calc_force!(ps, ϵ)  # Acceleration at one-step later
    
    for p in ps
        @. p.v += p.a * 0.5Δt  # Velocity at one-step later
    end
end

