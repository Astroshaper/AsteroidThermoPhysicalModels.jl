

function run_Hermite4(ps, params_sim)
    @unpack Δt, t_end, ϵ, α = params_sim
    
#     if haskey(params_sim, :Δt)
#         @unpack Δt = params_sim
#         Initialize!(ps, params_sim)
#     else
#         Initialize!(ps, ϵ)
#         Δt = get_Δt_Aarseth(⁺p, params_sim)
#     end

    times = (0:Δt:t_end)
    ⁺ps = deepcopy(ps)
    Initialize!(ps, ⁺ps, Δt, ϵ)

    ts, rs, vs, as = prep_snapshot(ps, params_sim)
    
    @show Δt
    @show t_end
    println(length(times), " timesteps (", length(ts), " to be saved)")
    
    snapshots = StructArray(Particle[])
        
    for (i, t) in enumerate(times)
        save_snapshot!(ts, rs, vs, as, i, t, ps, params_sim)
    
        predict!(ps, ⁺ps, Δt)
        evaluate!(ps, ⁺ps, ϵ)
        collect!(ps, ⁺ps, Δt, α)
        evaluate!(ps, ⁺ps, ϵ)

        prepare!(ps, ⁺ps, Δt)
    end
    ts, rs, vs, as
end



function Hermite4_shared()
    
    
end


#############################################################
#                      
#############################################################


function Initialize!(ps, ⁺ps, Δt, ϵ)
    evaluate!(ps, ⁺ps, ϵ)
    prepare!(ps, ⁺ps, Δt)
end


function predict!(ps, ⁺ps, Δt)
    Δt² = Δt * Δt
    Δt³ = Δt * Δt²

    for (p, ⁺p) in zip(ps, ⁺ps)
        r, v, a, a¹ = SVector{3}(p.r), SVector{3}(p.v), SVector{3}(p.a), SVector{3}(p.a¹)
        
        ⁺p.r .= r + Δt*v + Δt²/2*a + Δt³/6*a¹
        ⁺p.v .= v + Δt*a + Δt²/2*a¹
    end
end


function evaluate!(ps, ⁺ps, ϵ)
    for (pᵢ, ⁺pᵢ) in zip(ps, ⁺ps)
        ⁺pᵢ.a  .= 0.
        ⁺pᵢ.a¹ .= 0.
        for (pⱼ, ⁺pⱼ) in zip(ps, ⁺ps)
            pⱼ == pᵢ && continue
            pⱼ.m == 0 && continue 
            
            r = SVector{3}(⁺pᵢ.r) - SVector{3}(⁺pⱼ.r)
            v = SVector{3}(⁺pᵢ.v) - SVector{3}(⁺pⱼ.v)
            
            r² = norm(r)^2 + ϵ^2
            r⁻³ = r²^(-3/2)
            r⁻⁵ = r²^(-5/2)
            
            ⁺pᵢ.a  .-= G * pⱼ.m * r⁻³ * r
            ⁺pᵢ.a¹ .-= G * pⱼ.m * (r⁻³*v - 3*r⁻⁵*(v ⋅ r)*r)
        end
    end
end


function collect!(ps, ⁺ps, Δt, α)
    Δt² = Δt * Δt
    Δt³ = Δt * Δt²
    Δt⁴ = Δt * Δt³
    Δt⁵ = Δt * Δt⁴

    for (p, ⁺p) in zip(ps, ⁺ps)
         a,  a¹ = SVector{3}( p.a), SVector{3}( p.a¹)
        ⁺a, ⁺a¹ = SVector{3}(⁺p.a), SVector{3}(⁺p.a¹)
        
        p.a² .= ( -6(a - ⁺a) - Δt*(4a¹ + 2⁺a¹) ) / Δt²
        p.a³ .= ( 12(a - ⁺a) + Δt*6(a¹ +  ⁺a¹) ) / Δt³
        
        a², a³ = SVector{3}(p.a²), SVector{3}(p.a³)
                
        ⁺p.r .+= Δt⁴/24*a² + α*Δt⁵/120*a³
        ⁺p.v .+= Δt³/ 6*a² +   Δt⁴/ 24*a³
    end
end


function prepare!(ps, ⁺ps, Δt)
    for (p, ⁺p) in zip(ps, ⁺ps)
        @. p.r  = ⁺p.r
        @. p.v  = ⁺p.v
        @. p.a  = ⁺p.a
        @. p.a¹ = ⁺p.a¹
    
        # @. ⁺p.a³ = p.a³
        @. ⁺p.a² = p.a² + Δt*p.a³
    end
end


#############################################################
#                         Time step     
#############################################################


get_Δt_Aarseth!(⁺ps, params_sim) = minimum(get_Δt_Aarseth(⁺p, params_sim) for ⁺p in ⁺ps)


function get_Δt_Aarseth(⁺p::Particle, params_sim)
    @unpack η = params_sim

    s  = norm(⁺p.a)
    s¹ = norm(⁺p.a¹)
    s² = norm(⁺p.a²)
    s³ = norm(⁺p.a³)
    
    Δt = η * sqrt((s*s² + s¹*s¹) / (s¹*s³ + s²*s²))
end



