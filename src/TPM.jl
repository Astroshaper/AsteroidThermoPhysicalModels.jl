


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
# """
# function run_TPM!(shape::ShapeModel, orbit::OrbitalElements, spin::SpinParams, thermo_params::AbstractThermoParams, savepath="tmp.jld2")
#     @unpack P, Δt, t_begin, t_end = thermo_params
    
#     init_temps_zero!(shape, thermo_params)

#     ts = (t_begin:Δt:t_end) * P
#     timestamp = prep_timestamp(ts)
#     # surf_temp_table = zeros(length(shape.faces), Int(1/thermo_params.Δt)-1)

#     for (i, t) in enumerate(ts)
#         update_orbit!(orbit, t)
#         update_spin!(spin, t)
            
#         r̂☉ = normalize(orbit.r) * -1  # Shift the origin from the sun to the body
#         r̂☉ = orbit_to_body(r̂☉, spin)
        
#         update_flux_sun!(shape, orbit.F☉, r̂☉)
#         update_flux_scat_single!(shape, thermo_params)
#         update_flux_rad_single!(shape, thermo_params)
        
#         update_force!(shape, thermo_params)
#         sum_force_torque!(shape)
        
#         f = SVector{3}(shape.force)   # Body-fixed frame
#         τ = SVector{3}(shape.torque)  # Body-fixed frame

#         f = body_to_orbit(f, spin)  # Orbital plane frame
#         τ = body_to_orbit(τ, spin)  # Orbital plane frame

#         E_in, E_out, E_cons = energy_io(shape, thermo_params)

#         save_timestamp!(timestamp, i, orbit.u, orbit.ν, spin.ϕ, f..., τ..., E_in, E_out, E_cons)
        
#         update_temps!(shape, thermo_params)
#     end
#     mean_energy_cons_frac!(timestamp, spin)
#     jldsave(savepath; shape, orbit, spin, thermo_params, timestamp)

#     timestamp
# end

"""
- `shape`         : Shape model
- `et_range`      : Range of ephemeris times to run
- `sun`           : Sun's position in the body-fixed frame at epochs (Not normalized)
- `thermo_params` : Thermophysical parametes
- `savepath`      : Path to save data file
- `save_range`    : Indices in `et_range` to be saved

"""
function run_TPM!(shape::ShapeModel, et_range, sun, thermo_params::AbstractThermoParams, savepath, save_range)
 
    init_temps_zero!(shape, thermo_params)
    
    surf_temps = zeros(length(shape.faces), length(save_range))
    forces  = [zeros(3) for _ in eachindex(save_range)]
    torques = [zeros(3) for _ in eachindex(save_range)]
    
    idx = 1  # Index to save data

    for (et, r☉) in zip(et_range, sun)

        update_flux!(shape, r☉, thermo_params)

        if et_range[save_range[begin]] ≤ et ≤ et_range[save_range[end]]
            update_thermal_force!(shape, thermo_params)

            surf_temps[:, idx] .= surface_temperature(shape)
            forces[idx]  .= shape.force   # Body-fixed frame
            torques[idx] .= shape.torque  # Body-fixed frame
    
            idx += 1
        end

        E_in, E_out, E_cons = energy_io(shape, thermo_params)
        println(E_cons)

        update_temps!(shape, thermo_params)
    end
    
    jldsave(savepath; shape, et_range=et_range[save_range], sun=sun[save_range], thermo_params, surf_temps, forces, torques)
end

"""
    run_TPM!

Run TPM for a binary asteroid.

- shapes
- ephemerides
- thermo_params
- savepath
- savevalues
"""
function run_TPM!(shapes::Tuple, et_range, suns, S2P, d2_d1, thermo_params::AbstractThermoParams, savepath, savevalues)

    surf_temps = zeros(length(shapes[1].faces), length(et_range)), zeros(length(shapes[2].faces), length(et_range))
    forces  = [zeros(3) for _ in eachindex(et_range)], [zeros(3) for _ in eachindex(et_range)]
    torques = [zeros(3) for _ in eachindex(et_range)], [zeros(3) for _ in eachindex(et_range)]
    
    ## ProgressMeter setting
    # p = Progress(length(et_range); dt=1, desc="Running TPM...", showspeed=true)
    # ProgressMeter.ijulia_behavior(:clear)

    for (i, et) in enumerate(et_range)
        r☉₁ = suns[1][i]
        r☉₂ = suns[2][i]
        sec_from_pri = d2_d1[i]
        R₂₁ = S2P[i]

        ## Update enegey flux
        update_flux!(shapes[1], r☉₁, thermo_params)
        update_flux!(shapes[2], r☉₂, thermo_params)
        find_eclipse!(shapes, r☉₁, sec_from_pri, R₂₁)  # Mutual-shadowing

        ## Mutual-heating
        #
        #

        for (idx_shape, shape) in enumerate(shapes)
            update_thermal_force!(shape, thermo_params)

            surf_temps[idx_shape][:, i] .= surface_temperature(shape)
            forces[idx_shape][i]  .= shape.force   # Body-fixed frame
            torques[idx_shape][i] .= shape.torque  # Body-fixed frame
        end
    
        ## Energy input/output
        E_io_pri = energy_io(shapes[1], thermo_params)
        E_io_sec = energy_io(shapes[2], thermo_params)
        println(E_io_pri[3], ", ",  E_io_sec[3])

        ## Update the progress meter
        # showvalues = [(:i, i), (:E_cons_pri, E_io_pri[3]), (:E_cons_sec, E_io_sec[3])]
        # ProgressMeter.next!(p; showvalues)
        
        ## Update temperature distribution
        et == et_range[end] && break  # Stop to update the temperature at the final step
        update_temps!(shapes[1], thermo_params)
        update_temps!(shapes[2], thermo_params)
    end
    
    # jldsave(savepath; shapes, et_range, suns, S2P, thermo_params)
    jldsave(savepath; shapes, et_range, suns, S2P, thermo_params, surf_temps, forces, torques)
end


# ****************************************************************
#                     Convergence decision
# ****************************************************************

"""
    energy_io(shape::ShapeModel, params::AbstractThermoParams) -> E_in, E_out, E_cons

Input and output energy per second at a certain time

# Returns
- `E_in`   : Input energy per second at a certain time [W]
- `E_out`  : Output enegey per second at a certain time [W]
- `E_cons` : Output-input energy ratio (`E_out / E_in`)
"""
function energy_io(shape::ShapeModel, params::AbstractThermoParams)
    E_in   = energy_in(shape, params)
    E_out  = energy_out(shape, params)
    E_cons = E_out / E_in

    E_in, E_out, E_cons
end


"""
    energy_in(shape::ShapeModel, params::AbstractThermoParams) -> E_in
    energy_in(shape::ShapeModel, A_B)                          -> E_in
    energy_in(facet::Facet, A_B::Real, area)                   -> E_in

Input energy per second on a facet or a whole surface [W]
"""
energy_in(shape::ShapeModel, params::AbstractThermoParams) = energy_in(shape, params.A_B)

function energy_in(shape::ShapeModel, A_B)
    E_in = 0.
    for i in eachindex(shape.faces)
        E_in += energy_in(shape.facets[i], (A_B isa Real ? A_B : A_B[i]), shape.face_areas[i])
    end
    E_in
end

energy_in(facet::Facet, A_B::Real, area) = (1 - A_B) * (facet.flux.sun + facet.flux.scat) * area


"""
    energy_out(shape::ShapeModel, params::AbstractThermoParams) -> E_out
    energy_out(shape::ShapeModel, ε, A_TH)                      -> E_out
    energy_out(facet::Facet, ε::Real, A_TH::Real, area)         -> E_out

Output enegey per second from a facet or a whole surface [W]
"""
energy_out(shape::ShapeModel, params::AbstractThermoParams) = energy_out(shape, params.ε, params.A_TH)

function energy_out(shape::ShapeModel, ε, A_TH)
    E_out = 0.
    for i in eachindex(shape.faces)
        E_out += energy_out(shape.facets[i], (ε isa Real ? ε : ε[i]), (A_TH isa Real ? A_TH : A_TH[i]), shape.face_areas[i])
    end
    E_out
end

energy_out(facet::Facet, ε::Real, A_TH::Real, area) = ( ε*σ_SB*facet.temps[begin]^4 - (1 - A_TH)*facet.flux.rad ) * area


# ****************************************************************
#        Energy flux of sunlight, scattering, and radiation
# ****************************************************************


"""
    update_flux!(shape, r☉, thermo_params)

Update energy flux to every facet by solar radiation, scattering, and re-absorption of radiation
"""
function update_flux!(shape, r☉::AbstractVector, thermo_params)
    update_flux_sun!(shape, r☉)
    update_flux_scat_single!(shape, thermo_params)
    update_flux_rad_single!(shape, thermo_params)
end

"""
    update_flux_sun!(shape, r̂☉, F☉)

Update solar radiation flux on every facet of a shape model.

- `shape` : Shape model
- `r̂☉`    : Normalized vector indicating the direction of the sun in the body-fixed frame
- `F☉`    : Solar radiation flux [W/m²]
"""
function update_flux_sun!(shape, r̂☉, F☉)
    for i in eachindex(shape.faces)
        if isilluminated(shape, r̂☉, i)
            shape.facets[i].flux.sun = F☉ * (shape.face_normals[i] ⋅ r̂☉)
        else
            shape.facets[i].flux.sun = 0
        end
    end
end

"""
    update_flux_sun!(shape, r☉)

Update solar radiation flux on every facet of a shape model.

- `shape` : Shape model
- `r☉`    : Position of the sun in the body-fixed frame, which is not normalized.
"""
function update_flux_sun!(shape, r☉)
    r̂☉ = SVector{3}(normalize(r☉))
    F☉ = SOLAR_CONST / SPICE.convrt(norm(r☉), "m", "au")^2

    update_flux_sun!(shape, r̂☉, F☉)
end


"""
The secondary is within the critical angle to detect an eclipse event.
"""
function eclipse_is_possible(shapes, sun_from_pri, sec_from_pri)

    R₁ = maximum_radius(shapes[1])
    R₂ = maximum_radius(shapes[2])
    θ_crit = asin((R₁ + R₂) / norm(sec_from_pri))

    r̂☉ = SVector{3}(normalize(sun_from_pri))
    r̂ₛ = SVector{3}(normalize(sec_from_pri))
    θ = acos(r̂☉ ⋅ r̂ₛ)  # Angle of Sun-Primary-Secondary

    θ_crit < θ < π - θ_crit ? false : true
end


"""
    find_eclipse!(shapes, sun_from_pri, sec_from_pri, R₂₁)

- `shapes`
- `sun_from_pri`
- `sec_from_pri` : Position of the secondary relative to primary
- `R₂₁`          : Rotation matrix from secondary to primary
"""
function find_eclipse!(shapes, sun_from_pri, sec_from_pri, R₂₁)

    r̂☉ = SVector{3}(normalize(sun_from_pri))
    # r̂ₛ = SVector{3}(normalize(sec_from_pri))
    rₛ = SVector{3}(sec_from_pri)

    # cosθ = r̂☉ ⋅ r̂ₛ  # Cosine of angle of Sun-Primary-Secondary

    eclipse_is_possible(shapes, r̂☉, rₛ) == false && return

    for i in eachindex(shapes[1].faces)
        shapes[1].facets[i].flux.sun == 0 && continue
        A₁, B₁, C₁ = shapes[1].nodes[shapes[1].faces[i]]  # △ABC in primary
        G₁ = shapes[1].face_centers[i]                    # Center of △ABC in primary

        for j in eachindex(shapes[2].faces)
            shapes[2].facets[j].flux.sun == 0 && continue
            A₂, B₂, C₂ = shapes[2].nodes[shapes[2].faces[j]]  # △ABC in secondary
            G₂ = shapes[2].face_centers[j]                    # Center of △ABC in secondary
            
            ## Transform coordinates from secondary- to primary-fixed frame
            A₂ = R₂₁ * A₂ + rₛ
            B₂ = R₂₁ * B₂ + rₛ
            C₂ = R₂₁ * C₂ + rₛ
            G₂ = R₂₁ * G₂ + rₛ
            
            raycast(A₂, B₂, C₂, r̂☉, G₁) && (shapes[1].facets[i].flux.sun = 0)  # something wrong?
            raycast(A₁, B₁, C₁, r̂☉, G₂) && (shapes[2].facets[j].flux.sun = 0)  # something wrong?
        end
    end
end


"""
    update_flux_scat_single!(shape, params)
    update_flux_scat_single!(shape, A_B::Real)
    update_flux_scat_single!(shape, A_B::AbstractVector)

Single scattering of sunlight is considered.
"""
update_flux_scat_single!(shape, params::AbstractThermoParams) = update_flux_scat_single!(shape, params.A_B)

function update_flux_scat_single!(shape, A_B::Real)
    for facet in shape.facets
        facet.flux.scat = 0
        for visiblefacet in facet.visiblefacets
            facet.flux.scat += visiblefacet.f * A_B * shape.facets[visiblefacet.id].flux.sun
        end
    end
end

function update_flux_scat_single!(shape, A_B::AbstractVector)
    for facet in shape.facets
        facet.flux.scat = 0
        for visiblefacet in facet.visiblefacets
            facet.flux.scat += visiblefacet.f * A_B[visiblefacet.id] * shape.facets[visiblefacet.id].flux.sun
        end
    end
end

# """
#     update_flux_scat_mult!(shape, params)

# Multiple scattering of sunlight is considered.
# """
# update_flux_scat_mult!(shape, params::AbstractThermoParams) = update_flux_scat_mult!(shape, params.A_B)

# function update_flux_scat_mult!(shape, A_B::Real)
#     for facet in shape.facets
#         facet.flux.scat = 0
#         for (id, f) in zip(facet.visiblefacets.id, facet.visiblefacets.f)
#             facet.flux.scat += f * A_B * shape.facets[id].flux.sun
#         end
#     end
# end

# function update_flux_scat_mult!(shape, A_B::AbstractVector)
#     for facet in shape.facets
#         facet.flux.scat = 0
#         for (id, f) in zip(facet.visiblefacets.id, facet.visiblefacets.f)
#             facet.flux.scat += f * A_B[id] * shape.facets[id].flux.sun
#         end
#     end
# end

"""
    update_flux_rad_single!(shape, params)
    update_flux_rad_single!(shape, ε::Real)
    update_flux_rad_single!(shape, ε::AbstractVector)

Single radiation-reabsorption is considered,
assuming albedo is close to zero at thermal infrared wavelength.
"""
update_flux_rad_single!(shape, params::AbstractThermoParams) = update_flux_rad_single!(shape, params.ε)

function update_flux_rad_single!(shape, ε::Real)
    for facet in shape.facets
        facet.flux.rad = 0
        for visiblefacet in facet.visiblefacets
            T = shape.facets[visiblefacet.id].temps[begin]
            facet.flux.rad += visiblefacet.f * ε * σ_SB * T^4
        end
    end
end

function update_flux_rad_single!(shape, ε::AbstractVector)
    for facet in shape.facets
        facet.flux.rad = 0
        for visiblefacet in facet.visiblefacets
            T = shape.facets[visiblefacet.id].temps[begin]
            facet.flux.rad += visiblefacet.f * ε[visiblefacet.id] * σ_SB * T^4
        end
    end
end

