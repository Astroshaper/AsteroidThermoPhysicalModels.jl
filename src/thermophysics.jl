

# ****************************************************************
#              Thermal skin depth & Thermal inertia
# ****************************************************************

"""
    thermal_skin_depth(params)      -> l_2π
    thermal_skin_depth(P, k, ρ, Cp) -> l_2π

# Arguments
- `P`  : Cycle of thermal cycle [sec]
- `k`  : Thermal conductivity [W/m/K]
- `ρ`  : Material density [kg/m³]
- `Cp` : Heat capacity [J/kg/K]

# Return
- `l_2π` : Thermal skin depth [m]
"""
thermal_skin_depth(params) = thermal_skin_depth(params.P, params.k, params.ρ, params.Cp)
thermal_skin_depth(P, k, ρ, Cp) = @. √(4π * P * k / (ρ * Cp))


"""
    thermal_inertia(params)   -> Γ
    thermal_inertia(k, ρ, Cp) -> Γ

# Arguments
- `k`  : Thermal conductivity [W/m/K]
- `ρ`  : Material density [kg/m³]
- `Cp` : Heat capacity [J/kg/K]

# Return
- `Γ` : Thermal inertia [J ⋅ m⁻² ⋅ K⁻¹ ⋅ s⁻⁰⁵ (tiu)]
"""
thermal_inertia(params) = thermal_inertia(params.k, params.ρ, params.Cp)
thermal_inertia(k, ρ, Cp) = @. √(k * ρ * Cp)


# ****************************************************************
#               Struct for thermophysical properties
# ****************************************************************

"""
    struct ThermoParams{COMMON_INT, COMMON_FLOAT, FACET_INT, FACET_FLOAT}

# Fields
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

# Parametic types
- `COMMON_INT`   : Integer-typed property common for all facets
- `COMMON_FLOAT` : Float-typed property common for all facets
- `FACET_INT`    : Integer common for all facets or array giving an integer for each facet individually
- `FACET_FLOAT`  : Float common for all facets or array giving a float for each facet individually
"""
struct ThermoParams{COMMON_INT, COMMON_FLOAT, FACET_INT, FACET_FLOAT}
    A_B  ::FACET_FLOAT
    A_TH ::FACET_FLOAT
    k    ::FACET_FLOAT
    ρ    ::FACET_FLOAT
    Cp   ::FACET_FLOAT
    ϵ    ::FACET_FLOAT

    t_bgn::COMMON_FLOAT  # Common for all facets
    t_end::COMMON_FLOAT  # Common for all facets
    Δt   ::COMMON_FLOAT  # Common for all facets
    Nt   ::COMMON_INT    # Common for all facets

    z_max::FACET_FLOAT
    Δz   ::FACET_FLOAT
    Nz   ::FACET_INT

    P    ::COMMON_FLOAT  # Common for all facets
    l    ::FACET_FLOAT
    Γ    ::FACET_FLOAT
    λ    ::FACET_FLOAT
end


function ThermoParams(; A_B, A_TH, k, ρ, Cp, ϵ, t_bgn, t_end, Nt, z_max, Nz, P)

    t_bgn /= P                       # Normalized by period P
    t_end /= P                       # Normalized by period P
    Δt = (t_end - t_bgn) / (Nt - 1)  # Normalized by period P

    l = thermal_skin_depth(P, k, ρ, Cp)
    Γ = thermal_inertia(k, ρ, Cp)

    z_max = @. z_max / l      # Normalized by skin depth l
    Δz = @. z_max / (Nz - 1)  # Normalized by skin depth l

    λ = @. (Δt/Δz^2) / 4π
    maximum(λ) > 0.5 && println("λ should be smaller than 0.5 for convergence.")

    LENGTH = maximum(length.([A_B, A_TH, k, ρ, Cp, ϵ, z_max, Δz, Nz, l, Γ, λ]))

    if LENGTH > 1
        A_B   isa Real && (A_B   = fill(A_B,   LENGTH))
        A_TH  isa Real && (A_TH  = fill(A_TH,  LENGTH))
        k     isa Real && (k     = fill(k,     LENGTH))
        ρ     isa Real && (ρ     = fill(ρ,     LENGTH))
        Cp    isa Real && (Cp    = fill(Cp,    LENGTH))
        ϵ     isa Real && (ϵ     = fill(ϵ,     LENGTH))
        
        z_max isa Real && (z_max = fill(z_max, LENGTH))
        Δz    isa Real && (Δz    = fill(Δz,    LENGTH))
        Nz    isa Real && (Nz    = fill(Nz,    LENGTH))
        
        l     isa Real && (l     = fill(l,     LENGTH))
        Γ     isa Real && (Γ     = fill(Γ,     LENGTH))
        λ     isa Real && (λ     = fill(λ,     LENGTH))
    end
    
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
    update_temps!(shape::ShapeModel, params::ThermoParams)
    update_temps!(shape::ShapeModel, λ, A_B, A_TH, k, l, Δz, ϵ)
    update_temps!(facet::Facet, λ::Real, A_B::Real, A_TH::Real, k::Real, l::Real, Δz::Real, ϵ::Real)

Update temerature profie (`Facet.temps`) based on 1-D heat diffusion
"""
update_temps!(shape::ShapeModel, params::ThermoParams) = update_temps!(shape, params.λ, params.A_B, params.A_TH, params.k, params.l, params.Δz, params.ϵ)

function update_temps!(shape::ShapeModel, λ, A_B, A_TH, k, l, Δz, ϵ)
    step_heat_cond!(shape, λ)
    update_surf_temp!(shape, A_B, A_TH, k, l, Δz, ϵ)  # Surface boundary condition (Radiation)
    update_bottom_temp!(shape)                        # Internal boundary condition (Insulation)
end

function update_temps!(facet::Facet, λ::Real, A_B::Real, A_TH::Real, k::Real, l::Real, Δz::Real, ϵ::Real)
    step_heat_cond!(facet, λ)
    update_surf_temp!(facet, A_B, A_TH, k, l, Δz, ϵ)  # Surface boundary condition (Radiation)
    update_bottom_temp!(facet)                        # Internal boundary condition (Insulation)
end

"""
    step_heat_cond!(shape::ShapeModel, λ::AbstractVector)
    step_heat_cond!(shape::ShapeModel, λ::Real)
    step_heat_cond!(facet::Facet,      λ::Real)
    step_heat_cond!(Tⱼ::AbstractVector, Tⱼ₊₁::AbstractVector, λ::Real)

Calculate temperature profile at the next step and update `Facet.temps`

# Arguments
- `facet` : Surface facet (`Facet`)
- `λ`     : Coefficient of heat conduction equation
- `Tⱼ`     : Temperatures
- `Tⱼ₊₁`   : Temperatures at the next timestep
"""
function step_heat_cond!(shape::ShapeModel, λ::AbstractVector)
    for (i, facet) in enumerate(shape.facets)
        step_heat_cond!(facet, λ[i])
    end
end

function step_heat_cond!(shape::ShapeModel, λ::Real)
    for facet in shape.facets
        step_heat_cond!(facet, λ)
    end
end

step_heat_cond!(facet::Facet, λ::Real) = step_heat_cond!(facet.temps, facet._temps_, λ)

function step_heat_cond!(Tⱼ::AbstractVector, Tⱼ₊₁::AbstractVector, λ::Real)
    @. Tⱼ₊₁[begin+1:end-1] = @views (1-2λ)*Tⱼ[begin+1:end-1] + λ*(Tⱼ[begin+2:end] + Tⱼ[begin:end-2])
    @. Tⱼ = Tⱼ₊₁  # Update cells for next step
end


"""
    update_surf_temp!(shape::ShapeModel, A_B, A_TH, k, l, Δz, ϵ)
    update_surf_temp!(shape::ShapeModel, A_B::Real, A_TH::Real, k::Real, l::Real, Δz::Real, ϵ::Real)
    update_surf_temp!(facet::Facet,      A_B::Real, A_TH::Real, k::Real, l::Real, Δz::Real, ϵ::Real)
    update_surf_temp!(T::AbstractVector, F_total::Real, k::Real, l::Real, Δz::Real, ϵ::Real)

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
function update_surf_temp!(shape::ShapeModel, A_B, A_TH, k, l, Δz, ϵ)
    for (i, facet) in enumerate(shape.facets)
        update_surf_temp!(facet, A_B[i], A_TH[i], k[i], l[i], Δz[i], ϵ[i])
    end
end

function update_surf_temp!(shape::ShapeModel, A_B::Real, A_TH::Real, k::Real, l::Real, Δz::Real, ϵ::Real)
    for facet in shape.facets
        update_surf_temp!(facet, A_B, A_TH, k, l, Δz, ϵ)
    end
end

function update_surf_temp!(facet::Facet, A_B::Real, A_TH::Real, k::Real, l::Real, Δz::Real, ϵ::Real)
    F_total = flux_total(facet, A_B, A_TH)
    update_surf_temp!(facet.temps, F_total, k, l, Δz, ϵ)
end

function update_surf_temp!(T::AbstractVector, F_total::Real, k::Real, l::Real, Δz::Real, ϵ::Real)
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
Update bottom temperature under boundary condition of insulation
"""
function update_bottom_temp!(shape::ShapeModel)
    for facet in shape.facets
        update_bottom_temp!(facet)
    end
end

function update_bottom_temp!(facet::Facet)
    facet.temps[end] = facet.temps[end-1] 
end


"""
    flux_total(shape::Shape, A_B::AbstractVector, A_TH::AbstractVector) -> F_total
    flux_total(facet::Facet, A_B::Real,           A_TH::Real)           -> F_total

# Arguments
- `shape` : Shape model (`Shape`)
- `facet` : Surface facet (`Facet`)
- `A_B`   : Bond albedo
- `A_TH`  : Albedo in thermal infrared wavelength

Total energy absorbed by the facet
"""
flux_total(shape::ShapeModel, A_B::AbstractVector, A_TH::AbstractVector) = [flux_total(facet, A_B[i], A_TH[i]) for (i, facet) in enumerate(shape.facets)]

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

