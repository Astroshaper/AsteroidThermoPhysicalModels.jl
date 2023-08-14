


abstract type ThermoPhysicalModel end


"""
    struct SingleTPM <: ThermoPhysicalModel

#Fields
- `shape`          : Shape model
- `thermo_params`  : Thermophysical parameters

- `flux`           : Flux on each face. Matrix of size (Number of faces, 3). Three components are:
    - `flux[:, 1]`     : F_sun,  flux of direct sunlight
    - `flux[:, 2]`     : F_scat, flux of scattered light
    - `flux[:, 3]`     : F_rad,  flux of thermal emission from surrounding surface
- `temperature`    : 3D array in size of (Nz, Ns, Nt). Temperature according to depth cells (Nz), faces (Ns), and time steps in one periodic cycle (Nt).
    -         ⋅----------⋅
    -     Nt /          /|
    -       ⋅--- Ns ---⋅ |
    -       |          | |
    -    Nz |          | ⋅
    -       |          |/
    -       ⋅----------⋅

- `face_forces`    : Thermal force on each face
- `force`          : Thermal recoil force at body-fixed frame (Yarkovsky effect)
- `torque`         : Thermal recoil torque at body-fixed frame (YORP effect)

- `SELF_SHADOWING` : Flag to consider self-shadowing
- `SELF_HEATING`   : Flag to consider self-heating

# TO DO:
- 同時に time step と depth step に関するベクトルをもつのが良いかもしれない
- roughness_maps   ::ShapeModel[]
"""
struct SingleTPM <: ThermoPhysicalModel
    shape          ::ShapeModel
    thermo_params  ::Union{UniformThermoParams, NonUniformThermoParams}

    flux           ::Matrix{Float64}
    temperature    ::Array{Float64, 3}

    face_forces    ::Vector{SVector{3, Float64}}
    force          ::MVector{3, Float64}
    torque         ::MVector{3, Float64}

    SELF_SHADOWING ::Bool
    SELF_HEATING   ::Bool
end


function SingleTPM(shape, thermo_params, SELF_SHADOWING, SELF_HEATING)

    Nz = thermo_params.Nz
    Ns = length(shape.faces)
    Nt = thermo_params.Nt

    flux = zeros(Ns, 3)
    temperature = zeros(Nz, Ns, Nt)

    face_forces = [zero(SVector{3, Float64}) for _ in shape.faces]
    force  = zero(MVector{3, Float64})
    torque = zero(MVector{3, Float64})

    SingleTPM(shape, thermo_params, flux, temperature, face_forces, force, torque, SELF_SHADOWING, SELF_HEATING)
end


"""
    struct BinaryTPM <: ThermoPhysicalModel

#Fields
- `pri`              : TPM for the primary
- `sec`              : TPM for the secondary
- `MUTUAL_SHADOWING` : Flag to consider mutual shadowing
- `MUTUAL_HEATING`   : Flag to consider mutual heating
"""
struct BinaryTPM <: ThermoPhysicalModel
    pri              ::SingleTPM
    sec              ::SingleTPM

    MUTUAL_SHADOWING ::Bool
    MUTUAL_HEATING   ::Bool
end


# ****************************************************************
#                   Initialize temperatures
# ****************************************************************


"""
    subsolar_temperature(r☉) -> Tₛₛ

Subsolar temperature [K] on an asteroid at a heliocentric distance `r☉` [m],
assuming radiative equilibrium with zero conductivity.
"""
subsolar_temperature(r☉, params::AbstractThermoParams) = subsolar_temperature(r☉, params.A_B, params.ε)

function subsolar_temperature(r☉, A_B, ε)
    Φ = SOLAR_CONST / SPICE.convrt(norm(r☉), "m", "au")^2  # Energy flux at the solar distance [W/m²]
    Tₛₛ = ((1 - A_B) * Φ / (ε * σ_SB))^(1/4)

    return Tₛₛ
end


"""
    init_temperature_zero!(shape::ShapeModel, params::AbstractThermoParams)

Initialize all temperature cells at 0 K.
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


"""
    init_temperature!(shape::ShapeModel, params::AbstractThermoParams, T₀::Real)

Initialize all temperature cells at the given temperature `T₀`
"""
function init_temperature!(shape::ShapeModel, params::AbstractThermoParams, T₀::Real)
    init_temperature_zero!(shape, params)
    shape.temperature[:, :, :] .= T₀
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
function run_TPM!(shape::ShapeModel, thermo_params::AbstractThermoParams, ephem, savepath)
    
    surf_temps = zeros(length(shape.faces), length(ephem.time))
    forces  = [zeros(3) for _ in eachindex(ephem.time)]
    torques = [zeros(3) for _ in eachindex(ephem.time)]

    ## ProgressMeter setting
    p = Progress(length(ephem.time); dt=1, desc="Running TPM...", showspeed=true)
    ProgressMeter.ijulia_behavior(:clear)
    
    for nₜ in eachindex(ephem.time)
        et = ephem.time[nₜ]
        r☉ = ephem.sun[nₜ]

        update_flux!(shape, r☉, thermo_params, nₜ)
        update_thermal_force!(shape, thermo_params, nₜ)

        ## Save data
        surf_temps[:, nₜ] .= surface_temperature(shape, nₜ)
        forces[nₜ]  .= shape.force   # Body-fixed frame
        torques[nₜ] .= shape.torque  # Body-fixed frame

        E_in, E_out, E_cons = energy_io(shape, thermo_params, nₜ)

        ## Update the progress meter
        showvalues = [("Timestep", nₜ), ("E_cons = E_out / E_in", E_cons)]
        ProgressMeter.next!(p; showvalues)

        nₜ == length(ephem.time) && break  # Stop to update the temperature at the final step
        update_temperature!(shape, thermo_params, nₜ)
    end

    jldsave(savepath; shape, thermo_params, ephem, surf_temps, forces, torques)
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
function run_TPM!(shapes::Tuple, thermo_params::AbstractThermoParams, ephem, savepath)

    surf_temps = zeros(length(shapes[1].faces), length(ephem.time)), zeros(length(shapes[2].faces), length(ephem.time))
    forces  = [zeros(3) for _ in eachindex(ephem.time)], [zeros(3) for _ in eachindex(ephem.time)]
    torques = [zeros(3) for _ in eachindex(ephem.time)], [zeros(3) for _ in eachindex(ephem.time)]
    
    ## ProgressMeter setting
    p = Progress(length(ephem.time); dt=1, desc="Running TPM...", showspeed=true)
    ProgressMeter.ijulia_behavior(:clear)

    for nₜ in eachindex(ephem.time)
        r☉₁ = ephem.sun1[nₜ]
        r☉₂ = ephem.sun2[nₜ]
        sec_from_pri = ephem.sec[nₜ]
        R₂₁ = ephem.S2P[nₜ]

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
        E_cons_pri = energy_io(shapes[1], thermo_params, nₜ)[3]
        E_cons_sec = energy_io(shapes[2], thermo_params, nₜ)[3]

        ## Update the progress meter
        showvalues = [("Timestep", nₜ), ("E_cons for primary", E_cons_pri), ("E_cons for secondary", E_cons_sec)]
        ProgressMeter.next!(p; showvalues)
        
        ## Update temperature distribution
        nₜ == length(ephem.time) && break  # Stop to update the temperature at the final step
        update_temperature!(shapes[1], thermo_params, nₜ)
        update_temperature!(shapes[2], thermo_params, nₜ)
    end
    
    jldsave(savepath; shapes, thermo_params, ephem, surf_temps, forces, torques)
end

