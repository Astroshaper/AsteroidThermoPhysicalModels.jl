


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
    
    df = DataFrame()
    df[:, :t    ] = (t_bgn:Δt:t_end) * P
    df[:, :u    ] = zeros(Nt)
    df[:, :ϕ    ] = zeros(Nt)
    df[:, :f_x  ] = zeros(Nt)
    df[:, :f_y  ] = zeros(Nt)
    df[:, :f_z  ] = zeros(Nt)
    df[:, :τ_x  ] = zeros(Nt)
    df[:, :τ_y  ] = zeros(Nt)
    df[:, :τ_z  ] = zeros(Nt)
    df[:, :E_in ] = zeros(Nt)
    df[:, :E_out] = zeros(Nt)
    
    for (i, t) in enumerate(df.t)
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

        # df.t[i] = t
        df.u[i] = u
        df.ϕ[i] = ϕ

        df.f_x[i] = f[1]
        df.f_y[i] = f[2]
        df.f_z[i] = f[3]

        df.τ_x[i] = τ[1]
        df.τ_y[i] = τ[2]
        df.τ_z[i] = τ[3]

        df.E_in[i]  = energy_in(shape, params)
        df.E_out[i] = energy_out(shape, params)
        
        update_temps!(shape, params)
    end

    df[:, :ν] = u2ν(df.u, orbit)
    df[:, :E_cons] = df.E_out ./ df.E_in
    df[:, :Ē_cons] = [ mean(df.E_cons[@. row.t - spin.P ≤ df.t ≤ row.t]) for row in eachrow(df) ];
    
    df
end

# ****************************************************************
#                     Convergence decision
# ****************************************************************


"""
    energy_in(shape::Shape, params::ThermoParams) -> E_in
    energy_in(shape::Shape, A_B::AbstractVector ) -> E_in
    energy_in(shape::Shape, A_B::Real           ) -> E_in
    energy_in(facet::Facet, A_B::Real           ) -> E_in

Input energy per second at a certain time [W]
"""
energy_in(shape::Shape, params::ThermoParams) = energy_in(shape, params.A_B)
energy_in(shape::Shape, A_B::AbstractVector) = sum(energy_in(facet, A_B[i]) for (i, facet) in enumerate(shape.facets))
energy_in(shape::Shape, A_B::Real) = sum(energy_in(facet, A_B) for facet in shape.facets)
energy_in(facet::Facet, A_B::Real) = (1 - A_B) * (facet.flux.sun + facet.flux.scat) * facet.area


"""
    energy_out(shape::Shape, params::ThermoParams                   ) -> E_out
    energy_out(shape::Shape, ϵ::AbstractVector, A_TH::AbstractVector) -> E_out
    energy_out(shape::Shape, ϵ::Real,           A_TH::Real          ) -> E_out
    energy_out(facet::Facet, ϵ::Real,           A_TH::Real          ) -> E_out

Output enegey per second at a certain time [W]
"""
energy_out(shape::Shape, params::ThermoParams) = energy_out(shape, params.ϵ, params.A_TH)
energy_out(shape::Shape, ϵ::AbstractVector, A_TH::AbstractVector) = sum(energy_out(facet, ϵ[i], A_TH[i]) for (i, facet) in enumerate(shape.facets))
energy_out(shape::Shape, ϵ::Real, A_TH::Real) = sum(energy_out(facet, ϵ, A_TH) for facet in shape.facets)
energy_out(facet::Facet, ϵ::Real, A_TH::Real) = ( ϵ*σ_SB*facet.temps[begin]^4 - (1 - A_TH)*facet.flux.rad ) * facet.area

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
