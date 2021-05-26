

# c.f.
# https://github.com/nitadori/Hermite/blob/master/SRC/hermite4-k.cpp


#############################################################


function forward!(ps, Δt, ϵ, α, η)
    predict!(ps, Δt)
    evaluate_by_predictor!(ps, ϵ)
    collect!(ps, Δt, α)
    evaluate_by_corrector!(ps, ϵ)
    prepare!(ps, Δt, η)
end


function run_sim(ps, Δt, t_end, save_interval::Int, ϵ, α, η)
    
    nsteps = length(0:Δt:t_end)  # Number of time-steps for calculation
    
    times = collect(0:Δt:t_end)[begin:save_interval:end]  # Time-series to be saved
    nsaves = length(times)
    snapshots = Vector{typeof(ps)}(undef, nsaves)  # Snapshots of particles
    snapshots[begin] = deepcopy(ps)
            
    @show Δt
    @show t_end
    @show nsteps
    @show nsaves
    
    for i in 1:nsteps-1
        # t = Δt * (i - 1)
        forward!(ps, Δt, ϵ, α, η)
        i%save_interval == 0 && (snapshots[div(i, save_interval)+1] = deepcopy(ps))
    end
    times, snapshots
end


#############################################################
#                        Predict
#############################################################


function predict!(ps::AbstractVector, Δt)
    for p in ps
        predict!(p, Δt)
    end
end


function predict!(p::Particle, Δt)
    @. p.ᵖr = p.r + Δt*p.v  + (Δt^2)/2*p.a⁰ + (Δt^3)/6*p.a¹
    @. p.ᵖv = p.v + Δt*p.a⁰ + (Δt^2)/2*p.a¹
end


#############################################################
#                        Evaluate
#############################################################


function evaluate!(a⁰, a¹, r1, v1, r2, v2, m2, ϵ)
    r = r1 .- r2
    v = v1 .- v2
    
    r_norm = sqrt(norm(r)^2 + ϵ^2)  # Softening parameter ϵ
    r⁻³ = r_norm^-3
    r⁻⁵ = r_norm^-5
    
    @. a⁰ -= G * m2 * r⁻³ * r
    @. a¹ -= G * m2 * (r⁻³*v - 3*r⁻⁵*(v ⋅ r)*r)
end


evaluate_by_predictor!(p1, p2, ϵ) = evaluate!(p1.⁺a⁰, p1.⁺a¹, p1.ᵖr, p1.ᵖv, p2.ᵖr, p2.ᵖv, p2.m, ϵ)
evaluate_by_corrector!(p1, p2, ϵ) = evaluate!(p1.⁺a⁰, p1.⁺a¹, p1.ᶜr, p1.ᶜv, p2.ᶜr, p2.ᶜv, p2.m, ϵ)


function evaluate_by_predictor!(ps::AbstractVector, ϵ)
    for p in ps
        p.⁺a⁰ .= 0
        p.⁺a¹ .= 0
    end
    
    for p1 in ps
        for p2 in ps
            p1 != p2 && evaluate_by_predictor!(p1, p2, ϵ)
        end
    end
end


function evaluate_by_corrector!(ps::AbstractVector, ϵ)
    for p in ps
        p.⁺a⁰ .= 0
        p.⁺a¹ .= 0
    end
    
    for p1 in ps
        for p2 in ps
            p1 != p2 && evaluate_by_corrector!(p1, p2, ϵ)
        end
    end
end


function initialize!(ps::AbstractVector, ϵ)
    for p1 in ps
        p1.a⁰ .= 0
        p1.a¹ .= 0
        for p2 in ps
            p1 != p2 && evaluate!(p1.a⁰, p1.a¹, p1.r, p1.v, p2.r, p2.v, p2.m, ϵ)
        end
    end
end


#############################################################
#                        Collect
#############################################################


function collect!(ps::AbstractVector, Δt, α)
    for p in ps
        collect!(p, Δt, α)
    end
end


function collect!(p::Particle, Δt, α)
    @. p.a² = (-6*(p.a⁰ - p.⁺a⁰) - Δt*(4p.a¹ + 2p.⁺a¹)) / Δt^2
    @. p.a³ = (12*(p.a⁰ - p.⁺a⁰) + 6*Δt*(p.a¹ + p.⁺a¹)) / Δt^3

    @. p.ᶜr = p.ᵖr + (Δt^4)/24*p.a² + α*(Δt^5)/120*p.a³
    @. p.ᶜv = p.ᵖv + (Δt^3)/ 6*p.a² +   (Δt^4)/ 24*p.a³
end


#############################################################
#                 Prepare for the next step
#############################################################


prepare!(ps::AbstractVector, Δt, η) = minimum(prepare!(p, Δt, η) for p in ps)

function prepare!(p::Particle, Δt, η)
    p.r .= p.ᶜr
    p.v .= p.ᶜv
    
    p.a⁰ .= p.⁺a⁰
    p.a¹ .= p.⁺a¹
    
    ⁺a³ = p.a³
    ⁺a² = p.a² + Δt*p.a³
    
    ⁺Δt = get_Δt_Aarseth(p.⁺a⁰, p.⁺a¹, ⁺a², ⁺a³, η)
end


# get_Δt_Aarseth(p::Particle, η) = get_Δt_Aarseth(p.a⁰, p.a¹, p.a², p.a³, η)

function get_Δt_Aarseth(a⁰, a¹, a², a³, η)
    s⁰ = norm(a⁰)
    s¹ = norm(a¹)
    s² = norm(a²)
    s³ = norm(a³)
    
    η * sqrt((s⁰*s² + s¹*s¹) / (s¹*s³ + s²*s²))
end


get_Δt_initial(ps::AbstractVector, η) = minimum(get_Δt_initial(p, η) for p in ps)
get_Δt_initial(p::Particle, η) = get_Δt_initial(p.a⁰, p.a¹, η)
get_Δt_initial(a⁰, a¹, η) = η * norm(a⁰) / norm(a¹)





#############################################################

#############################################################


# function apply_force!(p, f)
#     p.a .+= f ./ p.m
# end


# function update!(ps, Δt)
#     for p in ps
#         update!(p, Δt)
#     end
# end


# function update!(p::Particle, Δt)
#     update_velocity!(p, Δt)
#     update_location!(p, Δt)
# end


# function update_velocity!(p, Δt)
#     p.v .+= p.a .* Δt
# end


# function update_location!(p, Δt)
#     p.r .+= p.v .* Δt
# end


# function reset_acceleration!(ps)
#     for p in ps
#         reset_acceleration!(p)
#     end
# end


# function reset_acceleration!(p::Particle)
#     p.a .= zero(eltype(p.a))
# end


###########################################

# function apply_mutual_gravity!(ps)
#     for p in ps
#         apply_net_gravity!(p, ps)
#     end
# end


# function apply_net_gravity!(p_i, ps)
#     for p_j in ps
#         p_i == p_j && continue
#         apply_gravity!(p_i, p_j)
#     end
# end


# function apply_gravity!(p1, p2)
#     Δr = get_distance(p1, p2)
#     f1 = G * p1.m * p2.m / Δr^3 .* (p2.r .- p1.r)
#     apply_force!(p1, f1)
# end


## 粒子の組み合わせを見つけて、ループで回すならこっちでもいいかも
# function get_gravity(p1, p2)
#     Δr = get_distance(p1, p2)
#     f1 = G * p1.m * p2.m / Δr^3 .* (p2.r .- p1.r)
#     f2 = f1 .* -1
    
#     f1, f2
# end

# get_distance(p1, p2) = norm(p1.r .- p2.r);