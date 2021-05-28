


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
    ϵ             = params_sim.ϵ
    Δt            = params_sim.Δt
    t_end         = params_sim.t_end
    save_interval = params_sim.save_interval
    
    times, snapshots = init_snapshot(ps, Δt, t_end, save_interval)
    
    for (i, t) in enumerate(times)
        save_snapshot!(i, save_interval, snapshots, ps)

        evaluate_Euler!(ps, ϵ)
        update!(ps, Δt)
    end

    times[begin:save_interval:end], snapshots
end


function init_snapshot(ps, Δt, t_end, save_interval)
    times = (0:Δt:t_end)
    snapshots = Vector{typeof(ps)}(undef, length(0:Δt*save_interval:t_end))
    times, snapshots
end
    

save_snapshot!(i, save_interval, snapshots, ps) = (i-1)%save_interval == 0 && (snapshots[i ÷ save_interval + 1] = deepcopy(ps))