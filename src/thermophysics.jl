

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

- `t_begin` : Start time of the thermophysical simulation [sec]
- `t_end`   : End time of the thermophysical simulation [sec]
- `Δt`      : Timestep [sec]
- `Nt`      : Number of timesteps

- `z_max` : Depth of the bottom of a heat conduction equation [m]
- `Δz`    : Depth step width [m]
- `Nz`    : Number of depth steps

- `P`     : Cycle of thermal cycle (rotation period) [sec]
- `l`     : Thermal skin depth [m]
- `Γ`     : Thermal inertia [J ⋅ m⁻² ⋅ K⁻¹ ⋅ s⁻⁰⁵ (tiu)]
- `λ`     : Non-dimensional coefficient for heat diffusion equation
"""
struct NonUniformThermoParams <: AbstractThermoParams
    A_B ::Vector{Float64}
    A_TH::Vector{Float64}
    k   ::Vector{Float64}
    ρ   ::Vector{Float64}
    Cp  ::Vector{Float64}
    ε   ::Vector{Float64}

    t_begin::Float64    # Common for all faces
    t_end  ::Float64    # Common for all faces
    Δt     ::Float64    # Common for all faces
    Nt     ::Int        # Common for all faces

    z_max::Float64      # Common for all faces
    Δz   ::Float64      # Common for all faces
    Nz   ::Int          # Common for all faces

    P::Float64          # Common for all faces
    l::Vector{Float64}
    Γ::Vector{Float64}
    λ::Vector{Float64}
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

- `t_begin` : Start time of the thermophysical simulation [sec]
- `t_end`   : End time of the thermophysical simulation [sec]
- `Δt`      : Timestep [sec]
- `Nt`      : Number of timesteps

- `z_max` : Depth of the bottom of a heat conduction equation [m]
- `Δz`    : Depth step width [m]
- `Nz`    : Number of depth steps

- `P`     : Thermal cycle (rotation period) [sec]
- `l`     : Thermal skin depth [m]
- `Γ`     : Thermal inertia [J ⋅ m⁻² ⋅ K⁻¹ ⋅ s⁻⁰⁵ (tiu)]
- `λ`     : Non-dimensional coefficient for heat diffusion equation
"""
struct UniformThermoParams <: AbstractThermoParams
    A_B ::Float64
    A_TH::Float64
    k   ::Float64
    ρ   ::Float64
    Cp  ::Float64
    ε   ::Float64

    t_begin::Float64
    t_end  ::Float64
    Δt     ::Float64
    Nt     ::Int

    z_max::Float64
    Δz   ::Float64
    Nz   ::Int

    P::Float64
    l::Float64
    Γ::Float64
    λ::Float64
end


"""
    thermoparams(; A_B, A_TH, k, ρ, Cp, ε, t_begin, t_end, Nt, z_max, Nz, P)
"""
function thermoparams(; A_B, A_TH, k, ρ, Cp, ε, t_begin, t_end, Nt, z_max, Nz, P)

    l = thermal_skin_depth(P, k, ρ, Cp)
    Γ = thermal_inertia(k, ρ, Cp)

    Δt = (t_end - t_begin) / (Nt - 1)
    Δz = z_max / (Nz - 1)

    λ = @. (Δt/P) / (Δz/l)^2 / 4π
    maximum(λ) > 0.5 && error("λ should be smaller than 0.5 for convergence of the forward Euler method.")

    LENGTH = maximum(length.([A_B, A_TH, k, ρ, Cp, ε, z_max, Δz, Nz, l, Γ, λ]))

    if LENGTH > 1
        A_B   isa Real && (A_B  = fill(A_B,  LENGTH))
        A_TH  isa Real && (A_TH = fill(A_TH, LENGTH))
        k     isa Real && (k    = fill(k,    LENGTH))
        ρ     isa Real && (ρ    = fill(ρ,    LENGTH))
        Cp    isa Real && (Cp   = fill(Cp,   LENGTH))
        ε     isa Real && (ε    = fill(ε,    LENGTH))
                
        l     isa Real && (l    = fill(l,    LENGTH))
        Γ     isa Real && (Γ    = fill(Γ,    LENGTH))
        λ     isa Real && (λ    = fill(λ,    LENGTH))
        
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
    
    msg =  "⋅-----------------------------------⋅\n"
    msg *= "|     Thermophysical parameters     |\n"
    msg *= "⋅-----------------------------------⋅\n"
    
    msg *= "  A_B     : $(A_B)\n"
    msg *= "  A_TH    : $(A_TH)\n"
    msg *= "  k       : $(k)\n"
    msg *= "  ρ       : $(ρ)\n"
    msg *= "  Cp      : $(Cp)\n"
    msg *= "  ε       : $(ε)\n"

    msg *= "-----------------------------------\n"
    msg *= "  t_begin : $(t_begin) [sec]\n"
    msg *= "            = $(t_begin / P) [P]\n"
    msg *= "  t_end   : $(t_end) [sec]\n"
    msg *= "            = $(t_end / P) [P]\n"
    msg *= "  Δt      : $(Δt) [sec]\n"
    msg *= "            = $(Δt / P) [P]\n"
    msg *= "  Nt      : $(Nt)\n"

    msg *= "-----------------------------------\n"
    msg *= "  z_max   : $(z_max) [m]\n"
    msg *= "            = $(z_max / l) [l]\n"
    msg *= "  Δz      : $(Δz) [m]\n"
    msg *= "            = $(Δz / l) [l]\n"
    msg *= "  Nz      : $(Nz)\n"
    
    msg *= "-----------------------------------\n"
    msg *= "  P       : $(P)\n"
    msg *= "  l       : $(l)\n"
    msg *= "  Γ       : $(Γ)\n"
    msg *= "  λ       : $(λ)\n"

    msg *= "-----------------------------------\n"
    print(io, msg)
end


# ****************************************************************
#                      1D heat conduction
# ****************************************************************

"""
    forward_temperature(shape::ShapeModel, λ, nₜ::Integer)

Calculate the temperature for the next time step (`nₜ + 1`) based on 1D heat conductivity equation.

TO DO: Allow selection of boundary conditions and solvers

# Arguments
- `shape`  : Shape model
- `params` : Thermophysical parameters
- `nₜ`     : Index of the current time step
"""
function update_temperature!(shape::ShapeModel, params::AbstractThermoParams, nₜ::Integer)
    λ = params.λ
    Tⱼ   = @views shape.temperature[:, :, nₜ  ]
    Tⱼ₊₁ = @views shape.temperature[:, :, nₜ+1]

    ## Forward Euler method
    @. Tⱼ₊₁[begin+1:end-1, :] = @views (1-2λ')*Tⱼ[begin+1:end-1, :] + λ'*(Tⱼ[begin+2:end, :] + Tⱼ[begin:end-2, :])

    ## Boundary conditions
    update_surface_temperature!(shape, params, nₜ+1)  # Radiation (surface)
    update_bottom_temperature!(shape, nₜ+1)           # Insulation (bottom)
end


# ****************************************************************
#                   Surface boundary condition
# ****************************************************************

"""
    update_surface_temperature!(shape::ShapeModel, params::AbstractThermoParams, nₜ::Integer)

Update surface temperature under radiative boundary condition using Newton's method

# Arguments
- `shape`  : Shape model (`ShapeModel`)
- `params` : Thermophysical prameters
- `nₜ`     : Index of the current time step

In the normalized equation of the surface boundary condition,
the coefficient `Γ / √(4π * P)` is equivalent for `k / l`,
where `Γ` is the thermal inertia and `P` the rotation period.
"""
function update_surface_temperature!(shape::ShapeModel, params::AbstractThermoParams, nₜ::Integer)
    for i in eachindex(shape.faces)
        F_sun, F_scat, F_rad = shape.flux[i, :]

        A_B  = (params.A_B  isa Real ? params.A_B  : params.A_B[i] )
        A_TH = (params.A_TH isa Real ? params.A_TH : params.A_TH[i])
        k    = (params.k    isa Real ? params.k    : params.k[i]   )
        l    = (params.l    isa Real ? params.l    : params.l[i]   )
        Δz   = (params.Δz   isa Real ? params.Δz   : params.Δz[i]  )
        ε    = (params.ε    isa Real ? params.ε    : params.ε[i]   )

        F_total = total_flux(A_B, A_TH, F_sun, F_scat, F_rad)
        update_surface_temperature!((@views shape.temperature[:, i, nₜ]), F_total, k, l, Δz/l, ε)  # Δz should be normalized by l
    end
end


"""
    update_surface_temperature!(T::AbstractVector, F_total::Real, k::Real, l::Real, Δz::Real, ε::Real)

Newton's method to update the surface temperature under radiative boundary condition

# Arguments
- `T`       : 1-D array of temperatures
- `F_total` : Total energy absorbed by the facet
- `k`       : Thermal conductivity [W/m/K]
- `l`       : Thermal skin depth [m]
- `Δz̄`      : Non-dimensional step in depth, normalized by thermal skin depth `l`
- `ε`       : Emissivity
"""
function update_surface_temperature!(T::AbstractVector, F_total::Real, k::Real, l::Real, Δz̄::Real, ε::Real)
    εσ = ε * σ_SB
    for _ in 1:20
        T_pri = T[begin]

        f = F_total + k / l * (T[begin+1] - T[begin]) / Δz̄ - εσ*T[begin]^4
        df = - k / l / Δz̄ - 4*εσ*T[begin]^3             
        T[begin] -= f / df

        err = abs(1 - T_pri / T[begin])
        err < 1e-10 && return
    end
end


# ****************************************************************
#                   Bottom boundary condition
# ****************************************************************

"""
    update_bottom_temperature!(shape::ShapeModel, nₜ::Integer)

Update bottom temperature under boundary condition of insulation
"""
function update_bottom_temperature!(shape::ShapeModel, nₜ::Integer)
    for i in eachindex(shape.faces)
        shape.temperature[end, i, nₜ] = shape.temperature[end-1, i, nₜ]
    end
end

