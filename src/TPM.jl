


# ****************************************************************
#                   Initialize temperatures
# ****************************************************************

"""
    init_temps_zero!(shape::ShapeModel, params)
    init_temps_zero!(shape::ShapeModel, Nz::AbstractVector)
    init_temps_zero!(shape::ShapeModel, Nz::Integer)
    init_temps_zero!(facet::Facet,      Nz::Integer)

Initialize temperature profile in depth on every facet.
All elements are intialized as 0 K.
"""
init_temps_zero!(shape::ShapeModel, params) = init_temps_zero!(shape, params.Nz)

function init_temps_zero!(shape::ShapeModel, Nz::AbstractVector)
    for (i, facet) in enumerate(shape.facets)
        init_temps_zero!(facet, Nz[i])
    end
end

function init_temps_zero!(shape::ShapeModel, Nz::Integer)
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

# """
#     mutable struct ThermoPhysicalModel

# # Fields
# - `shape`         :
# - `orbit`         :
# - `spin`          :
# - `thermo_params` : 

# - `t` : Time
# - `u` : Eccentric anomaly
# - `ν` : True anomaly
# - `ϕ` : Spin phase

# - `r`  : Position of the asteroid in the orbital plane frame
# - `F☉` : Solar irradiation at the position `r` [W/m²]
# - `r̂☉` : Unit vector for the solar position in the body-fixed frame
# """
# mutable struct ThermoPhysicalModel
#     shape        ::ShapeModel
#     orbit        ::OrbitalElements
#     spin         ::SpinParams
#     thermo_params::ThermoParams

#     t::Float64
#     u::Float64
#     ν::Float64
#     ϕ::Float64

#     r ::Vector{Float64}
#     F☉::Float64
#     r̂☉::Vector{Float64}
# end

# function ThermoPhysicalModel(shape, orbit, spin, thermo_params)
#     @unpack t_bgn, P = thermo_params

#     init_temps_zero!(shape, thermo_params)

#     t = t_bgn * P
#     u = solveKeplerEquation2(orbit, t)
#     ν = u2ν(u, orbit)
#     ϕ = spin.ω * t

#     r = get_r(orbit, u)
#     F☉ = getSolarIrradiation(norm(r))
    
#     r̂☉ = normalize(r) * -1  # Shift the origin from the sun to the body
#     r̂☉ = orbit_to_body(r̂☉, spin.γ, spin.ε, ϕ)

#     ThermoPhysicalModel(shape, orbit, spin, thermo_params, t, u, ν, ϕ, r, F☉, r̂☉)
# end

# function Base.show(io::IO, TPM::ThermoPhysicalModel)
#     @unpack t, u, ν, ϕ, r, F☉, r̂☉ = TPM
    
#     println("Time              : ", t)
#     println("Eccentric anomaly : ", rad2deg(u))
#     println("True anomaly      : ", rad2deg(ν))
#     println("Spin phase        : ", rad2deg(ϕ))

#     println("Body's position   : ", r)
#     println("Solar irradiation : ", F☉)
#     println("Sun's direction   : ", r̂☉)
# end

# """
#     nextstep!(TPM::ThermoPhysicalModel, Δt)

# Move on to the next timestep
# """
# function nextstep!(TPM::ThermoPhysicalModel, Δt)
#     @unpack t, u, ν, ϕ, r, F☉, r̂☉ = TPM

#     t += Δt
#     u = solveKeplerEquation2(orbit, t)
#     ν = u2ν(u, orbit)
#     ϕ = spin.ω * t

#     r = get_r(orbit, u)
#     F☉ = getSolarIrradiation(norm(r))
    
#     r̂☉ = normalize(r) * -1  # Shift the origin from the sun to the body
#     r̂☉ = orbit_to_body(r̂☉, spin.γ, spin.ε, ϕ)
# end

# """
# """
# function run(TPM::ThermoPhysicalModel)
#     @unpack shape, orbit, spin, thermo_params = 
#     @unpack P, Δt, t_bgn, t_end, Nt = params
    
    
#     init_temps_zero!(shape, params)

#     # ts = (t_bgn:Δt:t_end) * P
#     # outputs = [:u, :ν, :ϕ, :f_x, :f_y, :f_z, :τ_x, :τ_y, :τ_z, :E_in, :E_out, :E_cons]
    
#     for (i, t) in enumerate(ts)
#         @unapck shape = TPM

        
#         update_flux_sun!(shape, F☉, r̂☉)
#         update_flux_scat_single!(shape, params)
#         update_flux_rad_single!(shape, params)
        
#         update_force!(shape, params)
#         sum_force_torque!(shape)
        
#         f = SVector{3}(shape.force)   # Body-fixed frame
#         τ = SVector{3}(shape.torque)  # Body-fixed frame

#         f = body_to_orbit(f, spin.γ, spin.ε, ϕ)  # Orbital plane frame
#         τ = body_to_orbit(τ, spin.γ, spin.ε, ϕ)  # Orbital plane frame

#         E_in, E_out, E_cons = energy_io(shape, params)

#         values = [u, ν, ϕ, f..., τ..., E_in, E_out, E_cons]
#         save_to_dataframe!(df, i, outputs, values)
        
#         nextstep!(TPM, Δt)
#         update_temps!(shape, params)
#     end
#     # df[:, :Ē_cons] = [mean(df.E_cons[@. row.t - spin.P ≤ df.t ≤ row.t]) for row in eachrow(df)]
# end


"""
"""
function run_TPM(shape, orbit, spin, params::ThermoParams)
    @unpack P, Δt, t_bgn, t_end, Nt = params
    
    init_temps_zero!(shape, params)

    ts = (t_bgn:Δt:t_end) * P
    outputs = [:u, :ν, :ϕ, :f_x, :f_y, :f_z, :τ_x, :τ_y, :τ_z, :E_in, :E_out, :E_cons]
    df = prep_dataframe(ts, outputs)
    
    for (i, t) in enumerate(ts)
        ϕ = spin.ω * t

        u = solveKeplerEquation2(orbit, t)
        ν = u2ν(u, orbit)
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

        E_in, E_out, E_cons = energy_io(shape, params)

        values = [u, ν, ϕ, f..., τ..., E_in, E_out, E_cons]
        save_to_dataframe!(df, i, outputs, values)
        
        update_temps!(shape, params)
    end

    if "E_cons" in names(df)
        df[:, :Ē_cons] = [mean(df.E_cons[@. row.t - spin.P ≤ df.t ≤ row.t]) for row in eachrow(df)]
    end
    df
end

"""
"""
function prep_dataframe(ts, outputs; dtype=Float64)
    Nt = length(ts)
    df = DataFrame(t=ts)
    for symbol in outputs
        df[:, symbol] = Vector{dtype}(undef, Nt)
    end
    df
end

"""
"""
function save_to_dataframe!(df, i, outputs, values)
    for (symbol, value) in zip(outputs, values)
        df[i, symbol] = value
    end
end


# ****************************************************************
#                     Convergence decision
# ****************************************************************


"""
    energy_io(shape::ShapeModel, params::ThermoParams) -> E_in, E_out, E_cons

Input and output energy per second at a certain time

# Returns
- `E_in`   : Input energy per second at a certain time [W]
- `E_out`  : Output enegey per second at a certain time [W]
- `E_cons` : Output-input energy ratio (`E_out / E_in`)
"""
function energy_io(shape::ShapeModel, params::ThermoParams)
    E_in   = energy_in(shape, params)
    E_out  = energy_out(shape, params)
    E_cons = E_out / E_in

    E_in, E_out, E_cons
end


"""
    energy_in(shape::ShapeModel, params::ThermoParams) -> E_in
    energy_in(shape::ShapeModel, A_B::AbstractVector ) -> E_in
    energy_in(shape::ShapeModel, A_B::Real           ) -> E_in
    energy_in(facet::Facet,      A_B::Real           ) -> E_in

Input energy per second at a certain time [W]
"""
energy_in(shape::ShapeModel, params::ThermoParams) = energy_in(shape, params.A_B)
energy_in(shape::ShapeModel, A_B::AbstractVector) = sum(energy_in(facet, A_B[i]) for (i, facet) in enumerate(shape.facets))
energy_in(shape::ShapeModel, A_B::Real) = sum(energy_in(facet, A_B) for facet in shape.facets)
energy_in(facet::Facet, A_B::Real) = (1 - A_B) * (facet.flux.sun + facet.flux.scat) * facet.area


"""
    energy_out(shape::ShapeModel, params::ThermoParams                   ) -> E_out
    energy_out(shape::ShapeModel, ϵ::AbstractVector, A_TH::AbstractVector) -> E_out
    energy_out(shape::ShapeModel, ϵ::Real,           A_TH::Real          ) -> E_out
    energy_out(facet::Facet,      ϵ::Real,           A_TH::Real          ) -> E_out

Output enegey per second at a certain time [W]
"""
energy_out(shape::ShapeModel, params::ThermoParams) = energy_out(shape, params.ϵ, params.A_TH)
energy_out(shape::ShapeModel, ϵ::AbstractVector, A_TH::AbstractVector) = sum(energy_out(facet, ϵ[i], A_TH[i]) for (i, facet) in enumerate(shape.facets))
energy_out(shape::ShapeModel, ϵ::Real, A_TH::Real) = sum(energy_out(facet, ϵ, A_TH) for facet in shape.facets)
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
