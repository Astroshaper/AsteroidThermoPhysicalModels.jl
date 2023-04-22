


# ****************************************************************
#                      Analysis of YORP effect
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

