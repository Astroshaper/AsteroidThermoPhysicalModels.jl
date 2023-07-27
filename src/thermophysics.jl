

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

abstract type AbstractThermoParams end

"""
    struct NonUniformThermoParams

# Fields
- `A_B`   : Bond albedo
- `A_TH`  : Albedo at thermal radiation wavelength
- `k`     : Thermal conductivity [W/m/K]
- `ρ`     : Material density [kg/m³]
- `Cp`    : Heat capacity [J/kg/K]
- `ε`     : Emissivity

- `t_begin` : Start time of the simulation, normalized by period `P`
- `t_end`   : End time of the simulation, normalized by period `P`
- `Δt`      : Non-dimensional timesteps, normalized by period `P`
- `Nt`      : Number of timesteps

- `z_max` : Maximum depth for thermophysical simualtion, normalized by thermal skin depth `l`
- `Δz`    : Non-dimensional step in depth, normalized by thermal skin depth `l`
- `Nz`    : Number of depth steps

- `P`     : Cycle of thermal cycle (rotation period) [sec]
- `l`     : Thermal skin depth [m]
- `Γ`     : Thermal inertia [J ⋅ m⁻² ⋅ K⁻¹ ⋅ s⁻⁰⁵ (tiu)]
- `λ`     : Non-dimensional coefficient for heat diffusion equation
"""
struct NonUniformThermoParams <: AbstractThermoParams
    A_B  ::Vector{Float64}
    A_TH ::Vector{Float64}
    k    ::Vector{Float64}
    ρ    ::Vector{Float64}
    Cp   ::Vector{Float64}
    ε    ::Vector{Float64}

    t_begin::Float64  # Common for all facets
    t_end::Float64  # Common for all facets
    Δt   ::Float64  # Common for all facets
    Nt   ::Int    # Common for all facets

    z_max::Vector{Float64}
    Δz   ::Vector{Float64}
    Nz   ::Vector{Int}

    P    ::Float64  # Common for all facets
    l    ::Vector{Float64}
    Γ    ::Vector{Float64}
    λ    ::Vector{Float64}
end

"""
    struct UniformThermoParams

# Fields
- `A_B`   : Bond albedo
- `A_TH`  : Albedo at thermal radiation wavelength
- `k`     : Thermal conductivity [W/m/K]
- `ρ`     : Material density [kg/m³]
- `Cp`    : Heat capacity [J/kg/K]
- `ε`     : Emissivity

- `t_begin` : Start time of the simulation, normalized by period `P`
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
struct UniformThermoParams <: AbstractThermoParams
    A_B  ::Float64
    A_TH ::Float64
    k    ::Float64
    ρ    ::Float64
    Cp   ::Float64
    ε    ::Float64

    t_begin::Float64  # Common for all facets
    t_end::Float64  # Common for all facets
    Δt   ::Float64  # Common for all facets
    Nt   ::Int    # Common for all facets

    z_max::Float64
    Δz   ::Float64
    Nz   ::Int

    P    ::Float64  # Common for all facets
    l    ::Float64
    Γ    ::Float64
    λ    ::Float64
end


function thermoparams(; A_B, A_TH, k, ρ, Cp, ε, t_begin, t_end, Nt, z_max, Nz, P)
    t_begin /= P                       # Normalized by period P
    t_end /= P                       # Normalized by period P
    Δt = (t_end - t_begin) / (Nt - 1)  # Normalized by period P

    l = thermal_skin_depth(P, k, ρ, Cp)
    Γ = thermal_inertia(k, ρ, Cp)

    z_max = @. z_max / l      # Normalized by skin depth l
    Δz = @. z_max / (Nz - 1)  # Normalized by skin depth l

    λ = @. (Δt/Δz^2) / 4π
    maximum(λ) > 0.5 && println("λ should be smaller than 0.5 for convergence.")

    LENGTH = maximum(length.([A_B, A_TH, k, ρ, Cp, ε, z_max, Δz, Nz, l, Γ, λ]))

    if LENGTH > 1
        A_B   isa Real && (A_B   = fill(A_B,   LENGTH))
        A_TH  isa Real && (A_TH  = fill(A_TH,  LENGTH))
        k     isa Real && (k     = fill(k,     LENGTH))
        ρ     isa Real && (ρ     = fill(ρ,     LENGTH))
        Cp    isa Real && (Cp    = fill(Cp,    LENGTH))
        ε     isa Real && (ε     = fill(ε,     LENGTH))
        
        z_max isa Real && (z_max = fill(z_max, LENGTH))
        Δz    isa Real && (Δz    = fill(Δz,    LENGTH))
        Nz    isa Real && (Nz    = fill(Nz,    LENGTH))
        
        l     isa Real && (l     = fill(l,     LENGTH))
        Γ     isa Real && (Γ     = fill(Γ,     LENGTH))
        λ     isa Real && (λ     = fill(λ,     LENGTH))
        NonUniformThermoParams(A_B, A_TH, k, ρ, Cp, ε, t_begin, t_end, Δt, Nt, z_max, Δz, Nz, P, l, Γ, λ)
    else
        UniformThermoParams(A_B, A_TH, k, ρ, Cp, ε, t_begin, t_end, Δt, Nt, z_max, Δz, Nz, P, l, Γ, λ)
    end
end


function Base.show(io::IO, params::AbstractThermoParams)
    @unpack A_B, A_TH, k, ρ, Cp, ε = params
    @unpack t_begin, t_end, Δt, Nt = params
    @unpack z_max, Δz, Nz          = params
    @unpack P, l, Γ, λ             = params
    
    msg = "Thermophysical parameters\n"
    msg *= "-------------------------\n"
    
    msg *= "A_B   : $(A_B)\n"
    msg *= "A_TH  : $(A_TH)\n"
    msg *= "k     : $(k)\n"
    msg *= "ρ     : $(ρ)\n"
    msg *= "Cp    : $(Cp)\n"
    msg *= "ε     : $(ε)\n"

    msg *= "-------------------------\n"
    msg *= "t_begin : $(t_begin * P)\n"
    msg *= "t_begin : $(t_begin), (Normalized by period P)\n"
    msg *= "t_end   : $(t_end * P)\n"
    msg *= "t_end   : $(t_end), (Normalized by period P)\n"
    msg *= "Nt      : $(Nt)\n"
    msg *= "Δt      : $(Δt * P)\n"
    msg *= "Δt      : $(Δt), (Normalized by period P)\n"

    msg *= "-------------------------\n"
    msg *= "z_max : $(z_max * l)\n"
    msg *= "z_max : $(z_max), (Normalized by skin depth l)\n"
    msg *= "Nz    : $(Nz)\n"
    msg *= "Δz    : $(Δz * l)\n"
    msg *= "Δz    : $(Δz), (Normalized by skin depth l)\n"
    
    msg *= "-------------------------\n"
    msg *= "P     : $(P)\n"
    msg *= "l     : $(l)\n"
    msg *= "Γ     : $(Γ)\n"
    msg *= "λ     : $(λ)\n"

    msg *= "-------------------------\n"
    print(io, msg)
end


# ****************************************************************
#                      1D heat conduction
# ****************************************************************

"""
    update_temps!(shape::ShapeModel, params::AbstractThermoParams)
    update_temps!(shape::ShapeModel, λ, A_B, A_TH, k, l, Δz, ε)
    update_temps!(facet::Facet, λ::Real, A_B::Real, A_TH::Real, k::Real, l::Real, Δz::Real, ε::Real)

Update temerature profie (`Facet.temps`) based on 1-D heat diffusion
"""
update_temps!(shape::ShapeModel, params::AbstractThermoParams) = update_temps!(shape, params.λ, params.A_B, params.A_TH, params.k, params.l, params.Δz, params.ε)

function update_temps!(shape::ShapeModel, λ, A_B, A_TH, k, l, Δz, ε)
    step_heat_cond!(shape, λ)
    update_surf_temp!(shape, A_B, A_TH, k, l, Δz, ε)  # Surface boundary condition (Radiation)
    update_bottom_temp!(shape)                        # Internal boundary condition (Insulation)
end

function update_temps!(facet::Facet, λ::Real, A_B::Real, A_TH::Real, k::Real, l::Real, Δz::Real, ε::Real)
    step_heat_cond!(facet, λ)
    update_surf_temp!(facet, A_B, A_TH, k, l, Δz, ε)  # Surface boundary condition (Radiation)
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
    update_surf_temp!(shape::ShapeModel, params::AbstractThermoParams)
    update_surf_temp!(shape::ShapeModel, A_B, A_TH, k, l, Δz, ε)
    update_surf_temp!(T::AbstractVector, F_total::Real, k::Real, l::Real, Δz::Real, ε::Real)

Update surface temperature under radiative boundary condition using Newton's method

# Arguments
- `shape`   : Shape model (`ShapeModel`)
- `A_B`     : Bond albedo
- `A_TH`    : Albedo in thermal infrared wavelength
- `k`       : Thermal conductivity
- `l`       : Thermal skin depth
- `Δz`      : Step width in depth direction (normalized by thermal skin depth `l`)
- `ε`       : Emissivity

- `T`       : 1-D array of temperatures
- `F_total` : Total energy absorbed by the facet

In the normalized equation of the surface boundary condition,
the coefficient `Γ / √(4π * P)` is equivalent for `k / l`,
where `Γ` is the thermal inertia and `P` the rotation period.
"""
update_surf_temp!(shape::ShapeModel, params::AbstractThermoParams) = update_surf_temp!(shape, params.A_B, params.A_TH, params.k, params.l, params.Δz, params.ε)


function update_surf_temp!(shape::ShapeModel, A_B, A_TH, k, l, Δz, ε)
    for i in eachindex(shape.faces)
        A_B = (A_B isa Real ? A_B : A_B[i])
        A_TH = (A_TH isa Real ? A_TH : A_TH[i])
        F_sun, F_scat, F_rad = shape.flux[i, :]
        k = (k isa Real ? k : k[i])
        l = (l isa Real ? l : l[i])
        Δz = (Δz isa Real ? Δz : Δz[i])
        ε = (ε isa Real ? ε : ε[i])

        F_total = total_flux(A_B, A_TH, F_sun, F_scat, F_rad)
        update_surf_temp!(shape.facets[i].temps, F_total, k, l, Δz, ε)
    end
end


function update_surf_temp!(T::AbstractVector, F_total::Real, k::Real, l::Real, Δz::Real, ε::Real)
    εσ = ε * σ_SB
    for _ in 1:20
        T_pri = T[begin]

        f = F_total + k / l * (T[begin+1] - T[begin]) / Δz - εσ*T[begin]^4
        df = - k / l / Δz - 4*εσ*T[begin]^3             
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
    total_flux(A_B, A_TH, F_sun, F_scat, F_rad) -> F_total

Total energy absorbed by the face
"""
total_flux(A_B, A_TH, F_sun, F_scat, F_rad) = (1 - A_B) * (F_sun + F_scat) + (1 - A_TH) * F_rad
