

function run_Hermite4(ps, params_sim)
    ϵ             = params_sim.ϵ
    α             = params_sim.α
    η             = params_sim.η
    t_end         = params_sim.t_end
    save_interval = params_sim.save_interval
    
    if haskey(params_sim, :Δt)
        Δt = params_sim.Δt
        Initialize!(ps, ϵ)
    else
        Δt = Initialize!(ps, ϵ, η)
    end
    
    @show Δt

    times = (0:Δt:t_end)
    ⁺ps = deepcopy(ps)

    num_snapshot = length(0:Δt*save_interval:t_end)
    times_save = Vector{Float64}(undef, num_snapshot)
    snapshots = Vector{typeof(ps)}(undef, num_snapshot)
    
    for (i, t) in enumerate(times)
        save_snapshot!(i, save_interval, times_save, t, snapshots, ps)
        
        predict!(ps, ⁺ps, Δt)
        evaluate!(ps, ⁺ps, ϵ)
        collect!(ps, ⁺ps, Δt, α)
        evaluate!(ps, ⁺ps, ϵ)
        prepare!(ps, ⁺ps, Δt, η)
    end
    times_save, snapshots
end


function save_snapshot!(i, save_interval, times_save, t, snapshots, ps)
    if (i-1)%save_interval == 0
        i_save = i ÷ save_interval + 1
        times_save[i_save] = t
        snapshots[i_save] = deepcopy(ps)
    end
end


#############################################################
#                      
#############################################################

        
"""
    Initialize!(ps, ϵ)

# Fields
- `ps`
- `ϵ`
"""
function Initialize!(ps, ϵ)
    for pᵢ in ps
        pᵢ.a  .= 0.
        pᵢ.a¹ .= 0.
        for pⱼ in ps
            pᵢ == pⱼ && continue

            r = pᵢ.r - pⱼ.r
            v = pᵢ.v - pⱼ.v

            r_norm = sqrt(norm(r)^2 + ϵ^2)
            r⁻³ = r_norm^(-3)
            r⁻⁵ = r_norm^(-5)

            @. pᵢ.a  -= G * pⱼ.m * r⁻³ * r
            @. pᵢ.a¹ -= G * pⱼ.m * (r⁻³*v - 3*r⁻⁵*(v ⋅ r)*r)         
        end
    end 
end

"""
- `ps`
- `ϵ`
- `η`
"""
function Initialize!(ps, ϵ, η)
    Initialize!(ps, ϵ)
    Δt = minimum(η * norm(p.a) / norm(p.a¹) for p in ps)
end

"""
- `ps`
- `⁺ps`
- `Δt`
"""
function predict!(ps, ⁺ps, Δt)
    for (p, ⁺p) in zip(ps, ⁺ps)
        @. ⁺p.r = p.r + Δt*p.v + (Δt^2)/2*p.a + (Δt^3)/6*p.a¹
        @. ⁺p.v = p.v + Δt*p.a + (Δt^2)/2*p.a¹
    end
end

"""
- `ps`  :
- `⁺ps` :
- `ϵ`   : Softening parameter
"""
function evaluate!(ps, ⁺ps, ϵ)
    for (pᵢ, ⁺pᵢ) in zip(ps, ⁺ps)
        ⁺pᵢ.a  .= 0.
        ⁺pᵢ.a¹ .= 0.
        for (pⱼ, ⁺pⱼ) in zip(ps, ⁺ps)
            pᵢ == pⱼ && continue
            
            r = ⁺pᵢ.r - ⁺pⱼ.r
            v = ⁺pᵢ.v - ⁺pⱼ.v

            r_norm = sqrt(norm(r)^2 + ϵ^2)
            r⁻³ = r_norm^-3
            r⁻⁵ = r_norm^-5
            
            @. ⁺pᵢ.a  -= G * pⱼ.m * r⁻³ * r
            @. ⁺pᵢ.a¹ -= G * pⱼ.m * (r⁻³*v - 3*r⁻⁵*(v ⋅ r)*r)
        end
    end
end

"""
- `ps`
- `⁺ps`
- `Δt`
- `α`
"""
function collect!(ps, ⁺ps, Δt, α)
    for (p, ⁺p) in zip(ps, ⁺ps)
        @. p.a² = (-6*(p.a - ⁺p.a) - Δt*(4p.a¹ + 2⁺p.a¹)) / Δt^2
        @. p.a³ = (12*(p.a - ⁺p.a) + 6*Δt*(p.a¹ + ⁺p.a¹)) / Δt^3
        
        @. ⁺p.r += (Δt^4)/24*p.a² + α*(Δt^5)/120*p.a³
        @. ⁺p.v += (Δt^3)/ 6*p.a² +   (Δt^4)/ 24*p.a³
    end
end

"""
- `ps`
- `⁺ps`
- `Δt`
- `η`
"""
prepare!(ps, ⁺ps, Δt, η) = minimum(prepare!(p, ⁺p, Δt, η) for (p, ⁺p) in zip(ps, ⁺ps))

"""
"""
function prepare!(p::Particle, ⁺p::Particle, Δt, η)
    @. p.r  = ⁺p.r
    @. p.v  = ⁺p.v
    @. p.a  = ⁺p.a
    @. p.a¹ = ⁺p.a¹
    
    @. ⁺p.a³ = p.a³
    @. ⁺p.a² = p.a² + Δt*p.a³
    
    Δt = get_Δt_Aarseth(⁺p, η)
end

"""
"""
function get_Δt_Aarseth(⁺p, η)
    s  = norm(⁺p.a)
    s¹ = norm(⁺p.a¹)
    s² = norm(⁺p.a²)
    s³ = norm(⁺p.a³)
    
    Δt = η * sqrt((s*s² + s¹*s¹) / (s¹*s³ + s²*s²))
end



