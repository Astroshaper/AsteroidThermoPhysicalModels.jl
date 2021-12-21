
function run_Hermite4(ps, params_sim, filepath)
    @unpack Δt, t_end, ϵ, α, save_interval = params_sim

    times = (0:Δt:t_end)
    initialize!(ps)
    
    @show Δt
    @show t_end
    println(length(times))
    
    open(filepath, "w") do io
        
        E₀ = getTotalEnergy(ps)
        for (i, t) in enumerate(times)
            (i-1)%save_interval == 0 && save_snapshot_txt(io, t, ps)
                        
            predict!(ps, Δt)
            evaluate!(ps, ϵ)
            collect!(ps, Δt)
            evaluate!(ps, ϵ)

            prepare!(ps, Δt)
        end
        E₁ = getTotalEnergy(ps)
        ΔE = (E₁ - E₀)/E₀
        @show E₀
        @show E₁
        @show ΔE
    end
end

#############################################################

function run_Hermite4_test(ps, params_sim, filepath)
    open(filepath, "w") do io
    
        @unpack η, t_end, save_interval = params_sim

        initialize!(ps)
        Δt = η * minimum(norm(p.a)/norm(p.a¹) for p in ps)
        t = 0.
        i = 0
        
        E₀ = getTotalEnergy(ps)
        while t < t_end
            i%save_interval == 0 && save_snapshot_txt(io, t, ps)
            collide(ps) != (0, 0) && break
            
            Hermite4_shared!(ps, params_sim, Δt)
            Δt = get_Δt_Aarseth(ps, params_sim)
            t += Δt
            i += 1
        end
        E₁ = getTotalEnergy(ps)
        ΔE = (E₁ - E₀)/E₀
        @show E₀
        @show E₁
        @show ΔE
        println("Collision: ", collide(ps))
    end
end



function Hermite4_shared!(ps, params_sim, Δt)
    # @unpack ϵ, α = params_sim
    
    predict!(ps, Δt)
    evaluate!(ps)
    collect!(ps, Δt)
    evaluate!(ps)

    prepare!(ps, Δt)
end


#############################################################
#                   Predict & collect
#############################################################


function predict!(ps, Δt)
    Δt² = Δt * Δt
    Δt³ = Δt * Δt²
    
    for p in ps
        r, v, a, a¹ = SVector{3}(p.r), SVector{3}(p.v), SVector{3}(p.a), SVector{3}(p.a¹)
        
        p.⁺r .= r + Δt*v + Δt²/2*a + Δt³/6*a¹
        p.⁺v .= v + Δt*a + Δt²/2*a¹
    end
end


function evaluate!(ps, ϵ)
    for p in ps
        p.⁺a  .= 0
        p.⁺a¹ .= 0
    end
    
    for i in eachindex(ps)
        for j in eachindex(ps)
            i == j && continue
            ps[j].m == 0 && continue
            
            r = SVector{3}(ps[j].⁺r) - SVector{3}(ps[i].⁺r)
            v = SVector{3}(ps[j].⁺v) - SVector{3}(ps[i].⁺v)
            
            r² = norm(r)^2 + ϵ^2
            r⁻³ = r²^(-3/2)
            r⁻⁵ = r²^(-5/2)
            
            ps[i].⁺a  .+= G * ps[j].m * r * r⁻³
            ps[i].⁺a¹ .+= G * ps[j].m * (v*r⁻³ - 3(v ⋅ r)*r*r⁻⁵)     
        end
    end
end


evaluate!(ps) = evaluate!(ps, 0)


function collect!(ps, Δt, α)
    Δt² = Δt * Δt
    Δt³ = Δt * Δt²
    Δt⁴ = Δt * Δt³
    Δt⁵ = Δt * Δt⁴

    for p in ps
         a,  a¹ = SVector{3}(p.a),  SVector{3}(p.a¹)
        ⁺a, ⁺a¹ = SVector{3}(p.⁺a), SVector{3}(p.⁺a¹)
        
        p.a² .= ( -6(a - ⁺a) - Δt*(4a¹ + 2⁺a¹) ) / Δt²
        p.a³ .= ( 12(a - ⁺a) + 6Δt*(a¹ +  ⁺a¹) ) / Δt³
        
        a², a³ = SVector{3}(p.a²), SVector{3}(p.a³)
        
        p.⁺r .+= Δt⁴/24*a² + α*Δt⁵/120*a³
        p.⁺v .+= Δt³/ 6*a² +   Δt⁴/ 24*a³
    end
end


collect!(ps, Δt) = collect!(ps, Δt, 7/6)


function prepare!(ps, Δt)
    for p in ps
        @. p.r  = p.⁺r
        @. p.v  = p.⁺v
        @. p.a  = p.⁺a
        @. p.a¹ = p.⁺a¹
    
        # @. p.⁺a³ = p.a³
        @. p.a² += Δt*p.a³
    end
end


#############################################################
#                         Time step     
#############################################################


get_Δt_Aarseth(ps, params_sim) = minimum(get_Δt_Aarseth(p, params_sim) for p in ps)


function get_Δt_Aarseth(p::AbstractParticle, params_sim)
    @unpack η = params_sim

    a  = norm(SVector{3}(p.a ))
    a¹ = norm(SVector{3}(p.a¹))
    a² = norm(SVector{3}(p.a²))
    a³ = norm(SVector{3}(p.a³))
    
    Δt = η * sqrt((a*a² + a¹*a¹) / (a¹*a³ + a²*a²))
end



