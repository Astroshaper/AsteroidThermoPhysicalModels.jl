


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

    flux           ::Matrix{Float64}    # (Ns, 3)
    temperature    ::Array{Float64, 3}  # (Nz, Ns, Nt)

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


surface_temperature(stpm::SingleTPM, nₜ::Integer) = stpm.temperature[begin, :, nₜ]
surface_temperature(stpm::SingleTPM) = stpm.temperature[begin, :, end] 


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
    init_temperature!(stpm::SingleTPM, T₀::Real)

Initialize all temperature cells at the given temperature `T₀`
"""
function init_temperature!(stpm::SingleTPM, T₀::Real)
    stpm.temperature[:, :, :] .= T₀
end


"""
    init_temperature!(btpm::BinaryTPM, T₀::Real)

Initialize all temperature cells at the given temperature `T₀`
"""
function init_temperature!(btpm::BinaryTPM, T₀::Real)
    btpm.pri.temperature[:, :, :] .= T₀
    btpm.sec.temperature[:, :, :] .= T₀
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
    run_TPM!(stpm::SingleTPM, ephem, savepath)

Run TPM for a single asteroid.

# Arguments
- `stpm`     : Thermophysical model for a single asteroid
- `ephem`    : Ephemerides
    - `ephem.time` : Ephemeris times
    - `ephem.sun`  : Sun's position in the asteroid-fixed frame (Not normalized)
- `savepath` : Path to save data file
"""
function run_TPM!(stpm::SingleTPM, ephem, savepath)
    
    surf_temps = zeros(length(stpm.shape.faces), length(ephem.time))
    forces  = [zeros(3) for _ in eachindex(ephem.time)]
    torques = [zeros(3) for _ in eachindex(ephem.time)]

    ## ProgressMeter setting
    p = Progress(length(ephem.time); dt=1, desc="Running TPM...", showspeed=true)
    ProgressMeter.ijulia_behavior(:clear)
    
    for nₜ in eachindex(ephem.time)
        r☉ = ephem.sun[nₜ]

        update_flux_sun!(stpm, r☉)
        update_flux_scat_single!(stpm)
        update_flux_rad_single!(stpm, nₜ)
        
        update_thermal_force!(stpm, nₜ)

        ## Save data
        surf_temps[:, nₜ] .= surface_temperature(stpm, nₜ)
        forces[nₜ]  .= stpm.force   # Body-fixed frame
        torques[nₜ] .= stpm.torque  # Body-fixed frame

        ## Energy input/output
        E_in, E_out, E_cons = energy_io(stpm, nₜ)
        
        ## Update the progress meter
        showvalues = [("Timestep", nₜ), ("E_cons = E_out / E_in", E_cons)]
        ProgressMeter.next!(p; showvalues)

        nₜ == length(ephem.time) && break  # Stop to update the temperature at the final step
        update_temperature!(stpm, nₜ)
    end

    jldsave(savepath; stpm, ephem, surf_temps, forces, torques)
end

"""
    run_TPM!(btpm::BinaryTPM, ephem, savepath)

Run TPM for a binary asteroid.

# Arguments
- `stpm`     : Thermophysical model for a binary asteroid
- `ephem`    : Ephemerides
    - `time` : Ephemeris times
    - `sun1` : Sun's position in the primary's frame
    - `sun2` : Sun's position in the secondary's frame
    - `sec`  : Secondary's position in the primary's frame
    - `P2S`  : Rotation matrix from primary to secondary frames
    - `S2P`  : Rotation matrix from secondary to primary frames
- `savepath` : Path to save data file
"""
function run_TPM!(btpm::BinaryTPM, ephem, savepath)

    surf_temps = zeros(length(btpm.pri.shape.faces), length(ephem.time)), zeros(length(btpm.sec.shape.faces), length(ephem.time))
    forces  = [zeros(3) for _ in eachindex(ephem.time)], [zeros(3) for _ in eachindex(ephem.time)]
    torques = [zeros(3) for _ in eachindex(ephem.time)], [zeros(3) for _ in eachindex(ephem.time)]
    
    ## ProgressMeter setting
    p = Progress(length(ephem.time); dt=1, desc="Running TPM...", showspeed=true)
    ProgressMeter.ijulia_behavior(:clear)

    for nₜ in eachindex(ephem.time)
        r☉₁ = ephem.sun1[nₜ]  # Sun's position in the primary's frame
        r☉₂ = ephem.sun2[nₜ]  # Sun's position in the secondary's frame
        rₛ  = ephem.sec[nₜ]   # Secondary's position in the primary's frame
        R₂₁ = ephem.S2P[nₜ]   # Rotation matrix from secondary to primary frames

        ## Update enegey flux
        update_flux_sun!(btpm, r☉₁, r☉₂)
        update_flux_scat_single!(btpm)
        update_flux_rad_single!(btpm, nₜ)

        mutual_shadowing!(btpm, r☉₁, rₛ, R₂₁)  # Mutual-shadowing (eclipse)
        mutual_heating!(btpm, nₜ, rₛ, R₂₁)         # Mutual-heating

        update_thermal_force!(btpm, nₜ)

        ## Save data for primary
        surf_temps[1][:, nₜ] .= surface_temperature(btpm.pri, nₜ)
        forces[1][nₜ]  .= btpm.pri.force   # Body-fixed frame
        torques[1][nₜ] .= btpm.pri.torque  # Body-fixed frame

        ## Save data for secondary
        surf_temps[2][:, nₜ] .= surface_temperature(btpm.sec, nₜ)
        forces[2][nₜ]  .= btpm.sec.force   # Body-fixed frame
        torques[2][nₜ] .= btpm.sec.torque  # Body-fixed frame
    
        ## Energy input/output
        E_cons_pri = energy_io(btpm.pri, nₜ)[3]
        E_cons_sec = energy_io(btpm.sec, nₜ)[3]

        ## Update the progress meter
        showvalues = [("Timestep", nₜ), ("E_cons for primary", E_cons_pri), ("E_cons for secondary", E_cons_sec)]
        ProgressMeter.next!(p; showvalues)
        
        ## Update temperature distribution
        nₜ == length(ephem.time) && break  # Stop to update the temperature at the final step
        update_temperature!(btpm, nₜ)
    end
    
    jldsave(savepath; btpm, ephem, surf_temps, forces, torques)
end

