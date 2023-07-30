


# ****************************************************************
#                   Initialize temperatures
# ****************************************************************

"""
    init_temperature_zero!(shape::ShapeModel, params::AbstractThermoParams)

Initialize all cells in the temperature array at 0 K.
"""
function init_temperature_zero!(shape::ShapeModel, params::AbstractThermoParams)
    Nz = params.Nz
    Ns = length(shape.faces)
    Nt = params.Nt

    if size(shape.temperature) == (0, 0, 0)
        shape.temperature = zeros(Nz, Ns, Nt)
    elseif size(shape.temperature) == (Nz, Ns, Nt)
        shape.temperature .= 0.
    else
        error("ShapeModel.temperature has a wrong size.")
    end
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
    
    surf_temps = zeros(length(shape.faces), length(save_range))
    forces  = [zeros(3) for _ in eachindex(save_range)]
    torques = [zeros(3) for _ in eachindex(save_range)]
    
    idx = 1  # Index to save data

    for nₜ in eachindex(et_range)
        et = et_range[nₜ]
        r☉ = sun[nₜ]

        update_flux!(shape, r☉, thermo_params, nₜ)

        if et_range[save_range[begin]] ≤ et ≤ et_range[save_range[end]]
            update_thermal_force!(shape, thermo_params, nₜ)

            surf_temps[:, idx] .= surface_temperature(shape, nₜ)
            forces[idx]  .= shape.force   # Body-fixed frame
            torques[idx] .= shape.torque  # Body-fixed frame
    
            idx += 1
        end

        E_in, E_out, E_cons = energy_io(shape, thermo_params, nₜ)
        println(E_cons)

        nₜ == length(et_range) && break  # Stop to update the temperature at the final step
        update_temperature!(shape, thermo_params, nₜ)
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

    for nₜ in eachindex(et_range)
        et = et_range[nₜ]
        r☉₁ = suns[1][nₜ]
        r☉₂ = suns[2][nₜ]
        sec_from_pri = d2_d1[nₜ]
        R₂₁ = S2P[nₜ]

        ## Update enegey flux
        update_flux!(shapes[1], r☉₁, thermo_params, nₜ)
        update_flux!(shapes[2], r☉₂, thermo_params, nₜ)
        find_eclipse!(shapes, r☉₁, sec_from_pri, R₂₁)  # Mutual-shadowing

        ## Mutual-heating
        #
        #

        for (idx_shape, shape) in enumerate(shapes)
            update_thermal_force!(shape, thermo_params, nₜ)

            surf_temps[idx_shape][:, nₜ] .= surface_temperature(shape, nₜ)
            forces[idx_shape][nₜ]  .= shape.force   # Body-fixed frame
            torques[idx_shape][nₜ] .= shape.torque  # Body-fixed frame
        end
    
        ## Energy input/output
        E_io_pri = energy_io(shapes[1], thermo_params, nₜ)
        E_io_sec = energy_io(shapes[2], thermo_params, nₜ)
        println(E_io_pri[3], ", ",  E_io_sec[3])

        ## Update the progress meter
        # showvalues = [(:i, i), (:E_cons_pri, E_io_pri[3]), (:E_cons_sec, E_io_sec[3])]
        # ProgressMeter.next!(p; showvalues)
        
        ## Update temperature distribution
        nₜ == length(et_range) && break  # Stop to update the temperature at the final step
        update_temperature!(shapes[1], thermo_params, nₜ)
        update_temperature!(shapes[2], thermo_params, nₜ)
    end
    
    # jldsave(savepath; shapes, et_range, suns, S2P, thermo_params)
    jldsave(savepath; shapes, et_range, suns, S2P, thermo_params, surf_temps, forces, torques)
end


# ****************************************************************
#                     Energy input/output
# ****************************************************************

"""
    total_flux(A_B, A_TH, F_sun, F_scat, F_rad) -> F_total

Total energy absorbed by the face
"""
total_flux(A_B, A_TH, F_sun, F_scat, F_rad) = (1 - A_B) * (F_sun + F_scat) + (1 - A_TH) * F_rad


"""
    energy_io(shape::ShapeModel, params::AbstractThermoParams, nₜ::Integer) -> E_in, E_out, E_cons

Input and output energy per second at a certain time

# Returns
- `E_in`   : Input energy per second at a certain time [W]
- `E_out`  : Output enegey per second at a certain time [W]
- `E_cons` : Output-input energy ratio (`E_out / E_in`)
"""
function energy_io(shape::ShapeModel, params::AbstractThermoParams, nₜ::Integer)
    E_in   = energy_in(shape, params)
    E_out  = energy_out(shape, params, nₜ)
    E_cons = E_out / E_in

    E_in, E_out, E_cons
end


"""
    energy_in(shape::ShapeModel, params::AbstractThermoParams) -> E_in
    energy_in(shape::ShapeModel, A_B)                          -> E_in
    energy_in(A_B, F_sun, F_scat, area)                        -> E_in

Input energy per second on a facet or a whole surface [W]
"""
energy_in(shape::ShapeModel, params::AbstractThermoParams) = energy_in(shape, params.A_B)

function energy_in(shape::ShapeModel, A_B)
    E_in = 0.
    for i in eachindex(shape.faces)
        F_sun = shape.flux[i, 1]
        F_scat = shape.flux[i, 2]
        aᵢ = shape.face_areas[i]

        E_in += energy_in((A_B isa Real ? A_B : A_B[i]), F_sun, F_scat, aᵢ)
    end
    E_in
end

energy_in(A_B, F_sun, F_scat, area) = (1 - A_B) * (F_sun + F_scat) * area


"""
    energy_out(shape::ShapeModel, params::AbstractThermoParams, nₜ::Integer) -> E_out
    energy_out(shape, ε, A_TH, nₜ)                                           -> E_out
    energy_out(ε, T, A_TH, F_rad, area)                                      -> E_out

Output enegey per second from a facet or a whole surface [W]
"""
energy_out(shape::ShapeModel, params::AbstractThermoParams, nₜ::Integer) = energy_out(shape, params.ε, params.A_TH, nₜ)

function energy_out(shape, ε, A_TH, nₜ)
    E_out = 0.
    for i in eachindex(shape.faces)
        Tᵢ = shape.temperature[begin, i, nₜ]
        F_rad = shape.flux[i, 3]
        aᵢ = shape.face_areas[i]

        E_out += energy_out((ε isa Real ? ε : ε[i]), Tᵢ, (A_TH isa Real ? A_TH : A_TH[i]), F_rad, aᵢ)
    end
    E_out
end

energy_out(ε, T, A_TH, F_rad, area) = (ε * σ_SB * T^4 - (1 - A_TH) * F_rad) * area


# ****************************************************************
#                  Energy flux: Sunlight
# ****************************************************************

"""
    update_flux!(shape::ShapeModel, r☉::AbstractVector, params::AbstractThermoParams, nₜ::Integer)

Update energy flux to every facet by solar radiation, scattering, and re-absorption of radiation
"""
function update_flux!(shape::ShapeModel, r☉::AbstractVector, params::AbstractThermoParams, nₜ::Integer)
    update_flux_sun!(shape, r☉)
    update_flux_scat_single!(shape, params)
    update_flux_rad_single!(shape, params, nₜ)
end

"""
    update_flux_sun!(shape, r̂☉, F☉)

Update solar radiation flux on every facet of a shape model.

- `shape` : Shape model
- `r̂☉`    : Normalized vector indicating the direction of the sun in the body-fixed frame
- `F☉`    : Solar radiation flux [W/m²]
"""
function update_flux_sun!(shape::ShapeModel, r̂☉, F☉)
    for i in eachindex(shape.faces)
        if isilluminated(shape, r̂☉, i)
            n̂ = shape.face_normals[i]
            shape.flux[i, 1] = F☉ * (n̂ ⋅ r̂☉)
        else
            shape.flux[i, 1] = 0.
        end
    end
end

"""
    update_flux_sun!(shape, r☉)

Update solar radiation flux on every facet of a shape model.

- `shape` : Shape model
- `r☉`    : Position of the sun in the body-fixed frame, which is not normalized.
"""
function update_flux_sun!(shape::ShapeModel, r☉)
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
        shapes[1].flux[i, 1] == 0 && continue  # something wrong?
        A₁, B₁, C₁ = shapes[1].nodes[shapes[1].faces[i]]  # △ABC in primary
        G₁ = shapes[1].face_centers[i]                    # Center of △ABC in primary

        for j in eachindex(shapes[2].faces)
            shapes[2].flux[j, 1] == 0 && continue  # something wrong?
            A₂, B₂, C₂ = shapes[2].nodes[shapes[2].faces[j]]  # △ABC in secondary
            G₂ = shapes[2].face_centers[j]                    # Center of △ABC in secondary
            
            ## Transform coordinates from secondary- to primary-fixed frame
            A₂ = R₂₁ * A₂ + rₛ
            B₂ = R₂₁ * B₂ + rₛ
            C₂ = R₂₁ * C₂ + rₛ
            G₂ = R₂₁ * G₂ + rₛ
            
            raycast(A₂, B₂, C₂, r̂☉, G₁) && (shapes[1].flux[i, 1] = 0)  # something wrong?
            raycast(A₁, B₁, C₁, r̂☉, G₂) && (shapes[2].flux[j, 1] = 0)  # something wrong?
        end
    end
end


# ****************************************************************
#                 Energy flux: Scattering
# ****************************************************************

"""
    update_flux_scat_single!(shape, params::AbstractThermoParams)
    update_flux_scat_single!(shape, A_B)

Update flux of scattered sunlight, only considering single scattering.
"""
update_flux_scat_single!(shape::ShapeModel, params::AbstractThermoParams) = update_flux_scat_single!(shape, params.A_B)

function update_flux_scat_single!(shape, A_B)
    for i in eachindex(shape.faces)
        shape.flux[i, 2] = 0.
        for visiblefacet in shape.facets[i].visiblefacets
            j = visiblefacet.id
            fᵢⱼ = visiblefacet.f
            A_B = (A_B isa Real ? A_B : A_B[j])
            shape.flux[i, 2] += fᵢⱼ * A_B * shape.flux[j, 1]
        end
    end
end


##= TODO: Implement update_flux_scat_mult! =##

# """
#     update_flux_scat_mult!(shape, params::AbstractThermoParams)
#     update_flux_scat_mult!(shape, A_B)

# Update flux of scattered sunlight, considering multiple scattering.
# """


# ****************************************************************
#                Energy flux: Thermal radiation
# ****************************************************************

"""
    update_flux_rad_single!(shape, params::AbstractThermoParams, nₜ::Integer)
    update_flux_rad_single!(shape, ε, A_TH, nₜ)

Update flux of absorption of thermal radiation from surrounding surface.
Single radiation-absorption is only considered, assuming albedo is close to zero at thermal infrared wavelength.

# Arguments
- `shape`  : Shape model
- `params` : Thermophysical parameters
- `nₜ`     : Index of time step
"""
update_flux_rad_single!(shape::ShapeModel, params::AbstractThermoParams, nₜ::Integer) = update_flux_rad_single!(shape, params.ε, params.A_TH, nₜ)

function update_flux_rad_single!(shape, ε, A_TH, nₜ)
    for i in eachindex(shape.faces)
        shape.flux[i, 3] = 0.
        for visiblefacet in shape.facets[i].visiblefacets
            j = visiblefacet.id
            fᵢⱼ = visiblefacet.f
            ε = (ε isa Real ? ε : ε[j])
            A_TH = (A_TH isa Real ? A_TH : A_TH[j])
            Tⱼ = shape.temperature[begin, j, nₜ]
            
            shape.flux[i, 3] += ε * σ_SB * (1 - A_TH) * fᵢⱼ * Tⱼ^4
        end
    end
end

