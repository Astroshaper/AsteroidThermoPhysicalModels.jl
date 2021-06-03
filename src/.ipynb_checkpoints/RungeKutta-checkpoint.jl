

function run_RungeKutta(ps, params_sim)
    ϵ             = params_sim.ϵ
    η             = params_sim.η
    Δt            = params_sim.Δt
    t_end         = params_sim.t_end
    save_interval = params_sim.save_interval


    times = (0:Δt:t_end)
    times_save, snapshots = prep_snapshot!(ps, Δt, t_end, save_interval)
    # df = DataFrame()
    
    @show Δt
    @show t_end
    println(length(times), " timesteps (", length(times_save), " to be saved)")
    
    for (i, t) in enumerate(times)
        eval_acceleration!(ps, ϵ)
        (i-1)%save_interval == 0 && save_snapshot!(i, save_interval, times_save, t, snapshots, ps)
        
        step_RungeKutta2!(ps, Δt)
    end
    times_save, snapshots
end


################################################################


function eval_acceleration!(pᵢ, pⱼ, ϵ)
    r = pⱼ.r .- pᵢ.r
    r⁻³ = (norm(r)^2 + ϵ^2)^(-3/2)
    
    @. pᵢ.a += G * pⱼ.m * r⁻³ * r
end


function eval_acceleration!(ps, ϵ)
    for pᵢ in ps
        pᵢ.a .= 0.
        for pⱼ in ps
            pᵢ != pⱼ && eval_acceleration!(pᵢ, pⱼ, ϵ)
        end
    end
end


f_dr(p) = p.v  # dr/dt = v
f_dv(p) = p.a  # dv/dt = a

f(y, t) = [y[2], ]

r(ys) = [ ys[2], - (k/m)*ys[1] - g ]


function step_RungeKutta1!(ps, Δt)
    for p in ps
        @. p.r += p.v * Δt
        @. p.v += p.a * Δt
    end
end


function step_RungeKutta2!(ps, Δt)
    for p in ps
        @. p.r += p.v * Δt + 0.5 * p.a * Δt^2
        @. p.v += p.a * Δt
    end
end


function update_RungeKutta4!(f, t, y, Δt)
    k₁ = Δt * f(t, y)
    k₂ = Δt * f(t + Δt/2, y + k₁/2)
    k₃ = Δt * f(t + Δt/2, y + k₂/2)
    k₄ = Δt * f(t + Δt, y + k₃)
    
    y += (k₁ + 2k₂ + 2k₃ + k₄) / 6
end


################################################################

y0 = [pi - 0.1; 0.0]

function pend(t, y, dy)
    dy[1] = y[2]
    dy[2] = (-b * y[2]) - (c * sin(y[1]))
end

function f_pend(y, t)
    return [y[2], (-b * y[2]) - (c * sin(y[1]))]
end