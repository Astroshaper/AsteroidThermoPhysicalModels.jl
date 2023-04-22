

# ****************************************************************
#                  Photon recoil force / torque
# ****************************************************************

"""
    update_force!(shape::ShapeModel, params)
    update_force!(shape::ShapeModel, A_B::Real,           ϵ::Real)
    update_force!(shape::ShapeModel, A_B::AbstractVector, ϵ::Real)
    update_force!(shape::ShapeModel, A_B::Real,           ϵ::AbstractVector)
    update_force!(shape::ShapeModel, A_B::AbstractVector, ϵ::AbstractVector)
    update_force!(facet::Facet,      A_B::Real,           ϵ::Real)

Update photon recoil force on every facet (dfᵢ)
"""
update_force!(shape::ShapeModel, params) = update_force!(shape, params.A_B, params.ϵ)

function update_force!(shape::ShapeModel, A_B::Real, ϵ::Real)
    for facet in shape.facets
        update_force!(facet, A_B, ϵ)
    end
end

function update_force!(shape::ShapeModel, A_B::AbstractVector, ϵ::Real)
    for (i, facet) in enumerate(shape.facets)
        update_force!(facet, A_B[i], ϵ)
    end
end

function update_force!(shape::ShapeModel, A_B::Real, ϵ::AbstractVector)
    for (i, facet) in enumerate(shape.facets)
        update_force!(facet, A_B, ϵ[i])
    end
end

function update_force!(shape::ShapeModel, A_B::AbstractVector, ϵ::AbstractVector)
    for (i, facet) in enumerate(shape.facets)
        update_force!(facet, A_B[i], ϵ[i])
    end
end

function update_force!(facet::Facet, A_B::Real, ϵ::Real)
    E = A_B * facet.flux.scat + ϵ * σ_SB * facet.temps[begin]^4

    @. facet.force = facet.normal
    for vf in facet.visiblefacets
        @. facet.force -= 3/2 * vf.f * vf.d̂
    end
    @. facet.force *= - 2/3 * E * facet.area / c₀
end

"""
    sum_force_torque!(shape::Shape)

Integrate the force and torque over the global surface
"""
function sum_force_torque!(shape::ShapeModel)
    shape.force  .= 0
    shape.torque .= 0
    for facet in shape.facets
        r  = SVector{3}(facet.center)
        r̂  = normalize(r)
        df = SVector{3}(facet.force)
        
        shape.force  .+= (r̂ ⋅ df) * r̂  # Photon recoil force
        shape.torque .+= r × df        # Photon recoil torque
    end
end

"""
    update_force_Rubincam!(shape, params_thermo)

Update photon recoil force on every facet (df) based on Rubincam (2000) approximation
"""
function update_force_Rubincam!(shape, params_thermo)
    shape.torque .= 0.
    for facet in shape.facets
        facet.force .= - 2 * facet.flux.sun * facet.area / (3*c₀) .* facet.normal
        shape.torque .+= facet.center × SVector{3}(facet.force)
    end
end


# ****************************************************************
#                            Analysis
# ****************************************************************

getτω(τ, spin) = τ ⋅ spin.ŝ
getτε(τ, spin) = τ ⋅ spin_perp_unit1(spin)
getτψ(τ, spin) = τ ⋅ spin_perp_unit2(spin)


"""
    analyze_YORP(df, spin, MOI)

# Parameters
- `df`   : Dataframe of timestamp
- `spin` : `SpinParams`
- `C`    : Moment of inertia
"""
function analyze_YORP(df, spin, MOI)
    τ̄ = [mean(df.τ_x), mean(df.τ_y), mean(df.τ_z)]

    τ̄_ω = τ̄ ⋅ spin.ŝ
    τ̄_ε = τ̄ ⋅ spin_perp_unit1(spin)
    τ̄_ψ = τ̄ ⋅ spin_perp_unit2(spin)

    ## [rad/sec/sec]
    ω̇  = τ̄_ω / MOI
    ωε̇ = τ̄_ε / MOI
    ωψ̇ = τ̄_ψ / MOI

    ## [deg/day/day]
    ω̇  = rad2deg(ω̇)  * (3600*24)^2
    ωε̇ = rad2deg(ωε̇) * (3600*24)^2
    ωψ̇ = rad2deg(ωψ̇) * (3600*24)^2

    (τ̄ = τ̄, τ̄_ω = τ̄_ω, τ̄_ε = τ̄_ε, τ̄_ψ = τ̄_ψ, ω̇ = ω̇, ωε̇ = ωε̇, ωψ̇ = ωψ̇)
end


"""
    YORP_timescale(P_start, P_end, ω̇) -> timescale

Caluculate the YORP time scale given a constant acceleration/deceleration

# Parameters
- `P_start` [hour]
- `P_end`   [hour]
- `ω̇`       [deg/day/day]
"""
function YORP_timescale(ω̇, P_start, P_end)
    P_start *= 3600  # [sec]
    P_end   *= 3600  # [sec]

    ω_start = 2π / P_start  # Initial spin rate [rad/sec]
    ω_end   = 2π / P_end    # Final spin rate [rad/sec]
    
    Δω = ω_end - ω_start    # [rad/sec]
    Δω = rad2deg(Δω)        # [deg/sec]
    Δω *= 3600*24           # [deg/day]
    
    timescale = Δω / ω̇      # [day]
    SPICE.convrt(timescale, "days", "years")  # [year]
end

