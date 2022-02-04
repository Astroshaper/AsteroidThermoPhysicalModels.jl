

# ****************************************************************
#                Thermal skin depth and inertia
# ****************************************************************

"""
    thermal_skin_depth(P, k, ρ, Cp) -> l_2π

# Arguments
- `P`  : Cycle of thermal cycle [sec]
- `k`  : Thermal conductivity [W/m/K]
- `ρ`  : Material density [kg/m³]
- `Cp` : Heat capacity [J/kg/K]

# Return
- `l_2π` : Thermal skin depth [m]
"""
thermal_skin_depth(P, k, ρ, Cp) = √(4π * P * k / (ρ * Cp))
thermal_skin_depth(params) = thermal_skin_depth(params.P, params.k, params.ρ, params.Cp)

"""
    thermal_inertia(k, ρ, Cp) -> Γ

# Arguments
- `k`  : Thermal conductivity [W/m/K]
- `ρ`  : Material density [kg/m³]
- `Cp` : Heat capacity [J/kg/K]

# Return
- `Γ` : Thermal inertia [J ⋅ m⁻² ⋅ K⁻¹ ⋅ s⁻⁰⁵ (tiu)]
"""
thermal_inertia(k, ρ, Cp) = √(k * ρ * Cp)
thermal_inertia(params) = thermal_inertia(params.k, params.ρ, params.Cp)


# ****************************************************************
#               Struct for thermophysical properties
# ****************************************************************

"""
    struct ThermoParams

- `A_B`   : Bond albedo
- `A_TH`  : Albedo at thermal radiation wavelength
- `k`     : Thermal conductivity [W/m/K]
- `ρ`     : Material density [kg/m³]
- `Cp`    : Heat capacity [J/kg/K]
- `ϵ`     : Emissivity

- `t_bgn` : Start time of the simulation, normalized by period `P`
- `t_end` : End time of the simulation, normalized by period `P`
- `Δt`    : Non-dimensional timesteps, normalized by period `P`
- `Nt`    : Number of timesteps

- `z_max` : Maximum depth for thermophysical simualtion, normalized by thermal skin depth `l`
- `Δz`    : Non-dimensional step in depth, normalized by thermal skin depth `l`
- `Nz`    : Number of depth steps

- `P`     : Cycle of thermal cycle (rotation period) [sec]
- `l`     : Thermal skin depth [m]
- `Γ`     : Thermal inertia [J ⋅ m⁻² ⋅ K⁻¹ ⋅ s⁻⁰⁵ (tiu)]
- `λ`     : Non-dimensional coefficient for heat diffusion equation
"""
struct ThermoParams{T1, T2, T3, T4, T5, T6, T7, T8, T9, T10, T11}
    A_B  ::T1
    A_TH ::T2
    k    ::T3
    ρ    ::T4
    Cp   ::T5
    ϵ    ::T6

    t_bgn::T7
    t_end::T7
    Δt   ::T7
    Nt   ::T8

    z_max::T9
    Δz   ::T9
    Nz   ::T10

    P    ::T7
    l    ::T11
    Γ    ::T11
    λ    ::T11
end

function ThermoParams(; A_B, A_TH, k, ρ, Cp, ϵ, t_bgn=0., t_end, Nt, z_max, Nz, P)

    t_bgn /= P                       # Normalized by period
    t_end /= P                       # Normalized by period
    Δt = (t_end - t_bgn) / (Nt - 1)  # Normalized by period

    l = thermal_skin_depth(P, k, ρ, Cp)
    Γ = thermal_inertia(k, ρ, Cp)

    z_max /= l             # Normalized by skin depth
    Δz = z_max / (Nz - 1)  # Normalized by skin depth

    λ = 1/4π * (Δt/Δz^2)
    λ > 0.5 && println("λ = ", λ, ", which should be smaller than 0.5 for convergence.")
    
    ThermoParams(A_B, A_TH, k, ρ, Cp, ϵ, t_bgn, t_end, Δt, Nt, z_max, Δz, Nz, P, l, Γ, λ) 
end


function Base.show(io::IO, params::ThermoParams)
    @unpack A_B, A_TH, k, ρ, Cp, ϵ = params
    @unpack t_bgn, t_end, Δt, Nt   = params
    @unpack z_max, Δz, Nz          = params
    @unpack P, l, Γ, λ             = params
    
    println("Thermophysical parameters")
    println("-------------------------")
    
    println("A_B   : ", A_B)
    println("A_TH  : ", A_TH)
    println("k     : ", k)
    println("ρ     : ", ρ)
    println("Cp    : ", Cp)
    println("ϵ     : ", ϵ)

    println("-------------------------")
    println("t_bgn : ", t_bgn * P)
    println("t_bgn : ", t_bgn, " (Normalized by period P)")
    println("t_end : ", t_end * P)
    println("t_end : ", t_end, " (Normalized by period P)")
    println("Nt    : ", Nt)
    println("Δt    : ", Δt * P)
    println("Δt    : ", Δt, " (Normalized by period P)")

    println("-------------------------")
    println("z_max : ", z_max * l)
    println("z_max : ", z_max, " (Normalized by skin depth l)")
    println("Nz    : ", Nz)
    println("Δz    : ", Δz * l)
    println("Δz    : ", Δz, " (Normalized by skin depth l)")
    
    println("-------------------------")
    println("P     : ", P)
    println("l     : ", l)
    println("Γ     : ", Γ)
    println("λ     : ", λ)

    println("-------------------------")
end


# ****************************************************************
#                      1D heat conduction
# ****************************************************************

"""
    update_temps!(shape, params)
    update_temps!(facet::Facet, λ, A_B, A_TH, Δz, Γ, P, ϵ)

Update temerature profie (`Facet.temps`) based on 1-D heat diffusion
"""
function update_temps!(shape, params)
    @unpack λ, A_B, A_TH, k, l, Δz, ϵ = params

    for facet in shape.facets
        update_temps!(facet, λ, A_B, A_TH, k, l, Δz, ϵ)
    end
end

function update_temps!(shape::Shape, λ, A_B, A_TH, k, l, Δz, ϵ)
    # update_temps!.(shape.facets, λ, A_B, A_TH, k, l, Δz, ϵ)
    for facet in shape.facets
        update_temps!(facet, λ, A_B, A_TH, k, l, Δz, ϵ)
    end
end

function update_temps!(facet::Facet, λ, A_B, A_TH, k, l, Δz, ϵ)
    step_heat_cond!(facet, λ)
    update_surf_temp!(facet, A_B, A_TH, k, l, Δz, ϵ)  # Surface boundary condition (Radiation)
    facet.temps[end] = facet.temps[end-1]             # Internal boundary condition (Insulation)
end

"""
    step_heat_cond!(facet::Facet, λ)
    step_heat_cond!(Tⱼ, Tⱼ₊₁, λ)

Calculate temperature profile at the next step and update `Facet.temps`

# Arguments
- `facet` : Surface facet (`Facet`)
- `λ`     : Coefficient of heat conduction equation
- `Tⱼ`     : Temperatures
- `Tⱼ₊₁`   : Temperatures at the next timestep
"""
step_heat_cond!(facet::Facet, λ) = step_heat_cond!(facet.temps, facet._temps_, λ)

function step_heat_cond!(Tⱼ, Tⱼ₊₁, λ)
    @. Tⱼ₊₁[begin+1:end-1] = @views (1-2λ)*Tⱼ[begin+1:end-1] + λ*(Tⱼ[begin+2:end] + Tⱼ[begin:end-2])
    @. Tⱼ = Tⱼ₊₁  # Update cells for next step
end

"""
    update_surf_temp!(facet::Facet, A_B, A_TH, k, l, Δz, ϵ)
    update_surf_temp!(T, F_total, k, l, Δz, ϵ)

Update surface temperature under radiative boundary condition using Newton's method

# Arguments
- `facet`   : surface facet (`Facet`)
- `A_B`     : Bond albedo
- `A_TH`    : Albedo in thermal infrared wavelength
- `k`       : Thermal conductivity
- `l`       : Thermal skin depth
- `Δz`      : Step width in depth direction (normalized by thermal skin depth `l`)
- `ϵ`       : Emissivity

- `T`       : 1-D array of temperatures
- `F_total` : Total energy absorbed by the facet

In the normalized equation of the surface boundary condition,
the coefficient `Γ / √(4π * P)` is equivalent for `k / l`,
where `Γ` is the thermal inertia and `P` the rotation period.
"""
function update_surf_temp!(facet::Facet, A_B, A_TH, k, l, Δz, ϵ)
    F_total = flux_total(facet, A_B, A_TH)
    update_surf_temp!(facet.temps, F_total, k, l, Δz, ϵ)
end

function update_surf_temp!(T, F_total, k, l, Δz, ϵ)
    ϵσ = ϵ * σ_SB
    for _ in 1:20
        T_pri = T[begin]

        f = F_total + k / l * (T[begin+1] - T[begin]) / Δz - ϵσ*T[begin]^4
        df = - k / l / Δz - 4*ϵσ*T[begin]^3             
        T[begin] -= f / df

        err = abs(1 - T_pri / T[begin])
        err < 1e-10 && return
    end
end

"""
    flux_total(facet::Facet, A_B::Real, A_TH::Real) -> F_total

# Arguments
- `facet` : surface facet (`Facet`)
- `A_B`   : Bond albedo
- `A_TH`  : Albedo in thermal infrared wavelength

Total energy absorbed by the facet
"""
function flux_total(facet::Facet, A_B::Real, A_TH::Real)
    F_sun  = facet.flux.sun
    F_scat = facet.flux.scat
    F_rad  = facet.flux.rad
    
    F_total = (1 - A_B)*(F_sun + F_scat) + (1 - A_TH)*F_rad
end


# ****************************************************************
#
# ****************************************************************


"""
    intensity(λ, T) -> I

Intensity of radiation at a wavelength λ and tempertature T
according to the Planck function
"""
function intensity(λ, T)
    h = 6.62607015e-34  # Planck constant [J⋅s]
    k = 1.380649e-23    # Boltzmann's constant [J/K]

    I = 2 * h * c^2 / λ^5 / (exp(h * c₀ / (λ * k * T)) - 1)
end


ν2λ(ν) = c₀ / ν
λ2ν(λ) = c₀ / λ

