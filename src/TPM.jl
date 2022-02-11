


# ****************************************************************
#                   Initialize temperatures
# ****************************************************************

"""
    init_temps_zero!(shape::Shape, params)
    init_temps_zero!(shape::Shape, Nz::AbstractVector)
    init_temps_zero!(shape::Shape, Nz::Integer)
    init_temps_zero!(facet::Facet, Nz::Integer)

Initialize temperature profile in depth on every facet.
All elements are intialized as 0 K.
"""
init_temps_zero!(shape::Shape, params) = init_temps_zero!(shape, params.Nz)

function init_temps_zero!(shape::Shape, Nz::AbstractVector)
    for (i, facet) in enumerate(shape.facets)
        init_temps_zero!(facet, Nz[i])
    end
end

function init_temps_zero!(shape::Shape, Nz::Integer)
    for facet in shape.facets
        init_temps_zero!(facet, Nz)
    end
end

function init_temps_zero!(facet::Facet, Nz::Integer)
    isempty(facet.temps)   ? append!(facet.temps,   zeros(Nz)) : facet.temps   .= 0
    isempty(facet._temps_) ? append!(facet._temps_, zeros(Nz)) : facet._temps_ .= 0
end


# ****************************************************************
#                    Thermophysical modeling
# ****************************************************************

"""
    mutable struct ThermoPhysicalModel

# Fields
- `shape`         :
- `time`          :
- `orbit`         :
- `true_anomaly`  :
- `spin`          :
- `spin_phase`    :
- `thermo_params` : 
"""
mutable struct ThermoPhysicalModel
    shape        ::Shape
    time         ::Float64
    orbit        ::OrbitalElements
    true_anomaly ::Float64
    spin         ::SpinParams
    spin_phase   ::Float64
    thermo_params::ThermoParams
end

"""
"""
function run_TPM(shape, orbit, spin, params::ThermoParams)
    @unpack P, Δt, t_bgn, t_end, Nt = params
    
    init_temps_zero!(shape, params)
    
    ts = (t_bgn:Δt:t_end) * P
    df = DataFrame(
        t   = zeros(Nt), u   = zeros(Nt), ϕ   = zeros(Nt),
        f_x = zeros(Nt), f_y = zeros(Nt), f_z = zeros(Nt), 
        τ_x = zeros(Nt), τ_y = zeros(Nt), τ_z = zeros(Nt), 
    )

    for (i, t) in enumerate(ts)
        ϕ = spin.ω * t

        u = solveKeplerEquation2(orbit, t)
        r = get_r(orbit, u)
        F☉ = getSolarIrradiation(norm(r))
    
        r̂☉ = normalize(r) * -1  # Shift the origin from the sun to the body
        r̂☉ = orbit_to_body(r̂☉, spin.γ, spin.ε, ϕ)
        
        update_flux_sun!(shape, F☉, r̂☉)
        update_flux_scat_single!(shape, params)
        update_flux_rad_single!(shape, params)
        
        update_force!(shape, params)
        sum_force_torque!(shape)
        
        f = SVector{3}(shape.force)   # Body-fixed frame
        τ = SVector{3}(shape.torque)  # Body-fixed frame

        f = body_to_orbit(f, spin.γ, spin.ε, ϕ)  # Orbital plane frame
        τ = body_to_orbit(τ, spin.γ, spin.ε, ϕ)  # Orbital plane frame

        df.t[i] = t
        df.u[i] = u
        df.ϕ[i] = ϕ

        df.f_x[i] = f[1]
        df.f_y[i] = f[2]
        df.f_z[i] = f[3]

        df.τ_x[i] = τ[1]
        df.τ_y[i] = τ[2]
        df.τ_z[i] = τ[3]
        
        update_temps!(shape, params)
    end
    df
end


# ****************************************************************
#        Energy flux of sunlight, scattering, and radiation
# ****************************************************************

"""
    update_flux_sun!(shape, F☉, r̂☉)

Update illumination.
"""
function update_flux_sun!(shape, F☉, r̂☉)
    for facet in shape.facets
        if isIlluminated(facet, r̂☉, shape)
            facet.flux.sun = F☉ * (facet.normal ⋅ r̂☉)
        else
            facet.flux.sun = 0
        end
    end
end

"""
    update_flux_scat_single!(shape, params)
    update_flux_scat_single!(shape, A_B::Real)
    update_flux_scat_single!(shape, A_B::AbstractVector)

Single scattering of sunlight is considered.
"""
update_flux_scat_single!(shape, params::ThermoParams) = update_flux_scat_single!(shape, params.A_B)

function update_flux_scat_single!(shape, A_B::Real)
    for facet in shape.facets
        facet.flux.scat = 0
        for (id, f) in zip(facet.visiblefacets.id, facet.visiblefacets.f)
            facet.flux.scat += f * A_B * shape.facets[id].flux.sun
        end
    end
end

function update_flux_scat_single!(shape, A_B::AbstractVector)
    for facet in shape.facets
        facet.flux.scat = 0
        for (id, f) in zip(facet.visiblefacets.id, facet.visiblefacets.f)
            facet.flux.scat += f * A_B[id] * shape.facets[id].flux.sun
        end
    end
end

"""
    update_flux_scat_mult!(shape, params)

Multiple scattering of sunlight is considered.
"""
update_flux_scat_mult!(shape, params::ThermoParams) = update_flux_scat_mult!(shape, params.A_B)

function update_flux_scat_mult!(shape, A_B::Real)
    # for facet in shape.facets
    #     facet.flux.scat = 0
    #     for (id, f) in zip(facet.visiblefacets.id, facet.visiblefacets.f)
    #         facet.flux.scat += f * A_B * shape.facets[id].flux.sun
    #     end
    # end
end

function update_flux_scat_mult!(shape, A_B::AbstractVector)
    # for facet in shape.facets
    #     facet.flux.scat = 0
    #     for (id, f) in zip(facet.visiblefacets.id, facet.visiblefacets.f)
    #         facet.flux.scat += f * A_B[id] * shape.facets[id].flux.sun
    #     end
    # end
end

"""
    update_flux_rad_single!(shape, params)
    update_flux_rad_single!(shape, ϵ::Real)
    update_flux_rad_single!(shape, ϵ::AbstractVector)

Single radiation-reabsorption is considered,
assuming albedo is close to zero at thermal infrared wavelength.
"""
update_flux_rad_single!(shape, params::ThermoParams) = update_flux_rad_single!(shape, params.ϵ)

function update_flux_rad_single!(shape, ϵ::Real)
    for facet in shape.facets
        facet.flux.rad = 0
        for (id, f) in zip(facet.visiblefacets.id, facet.visiblefacets.f)
            T = shape.facets[id].temps[begin]
            facet.flux.rad += f * ϵ * σ_SB * T^4
        end
    end
end

function update_flux_rad_single!(shape, ϵ::AbstractVector)
    for facet in shape.facets
        facet.flux.rad = 0
        for (id, f) in zip(facet.visiblefacets.id, facet.visiblefacets.f)
            T = shape.facets[id].temps[begin]
            facet.flux.rad += f * ϵ[id] * σ_SB * T^4
        end
    end
end
