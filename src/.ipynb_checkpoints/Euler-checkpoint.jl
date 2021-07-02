


"""

- `⁺pᵢ` : Particle i at next step
- `pᵢ`  : Particle i
- `pⱼ`  : Particle j
- `ϵ`   : Softening parameter
"""
function evaluate_Euler!(⁺aᵢ, pᵢ, pⱼ, ϵ)
    r = pⱼ.r .- pᵢ.r
    r⁻³ = (norm(r)^2 + ϵ^2)^(-3/2)
    
    @. ⁺aᵢ += G * pⱼ.m * r⁻³ * r
end


function evaluate_Euler!(ps, ϵ)
    for pᵢ in ps
        pᵢ.a .= 0.
        for pⱼ in ps
            pᵢ != pⱼ && evaluate_Euler!(pᵢ.a, pᵢ, pⱼ, ϵ)
        end
    end
end


function update!(p::Particle, Δt)
    @. p.v += p.a * Δt
    @. p.r += p.v * Δt
end


function update!(ps, Δt)
    for p in ps
        update!(p, Δt)
    end
end


function run_Euler(ps, params_sim)
    @unpack ϵ, Δt, t_end = params_sim
    
    times = (0:Δt:t_end)
    ts, rs, vs, as = prep_snapshot(ps, params_sim)
    
    for (i, t) in enumerate(times)
        evaluate_Euler!(ps, ϵ)
        save_snapshot!(ts, rs, vs, as, i, t, ps, params_sim)
        update!(ps, Δt)
    end
    ts, rs, vs, as
end

