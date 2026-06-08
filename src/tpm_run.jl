#=
tpm_run.jl

Temperature initialization and thermophysical model simulation functions.
=#

# ╔═══════════════════════════════════════════════════════════════════╗
# ║                   Initialize temperatures                         ║
# ╚═══════════════════════════════════════════════════════════════════╝


"""
    subsolar_temperature(r☉, R_vis, ε) -> Tₛₛ
    subsolar_temperature(r☉, params::AbstractThermoParams) -> Tₛₛ

Calculate the subsolar temperature on an asteroid at a given heliocentric distance.

# Arguments
- `r☉::Vector`                   : Sun's position vector in the asteroid's fixed frame [m]
- `R_vis::Real`                  : Visible light reflectance (albedo) [-]
- `ε::Real`                      : Emissivity [-]
- `params::AbstractThermoParams` : Thermal parameters (alternative input)

# Returns
- `Tₛₛ::Float64` : Subsolar point temperature [K]

# Mathematical Formula
Assuming radiative equilibrium with zero thermal conductivity (zero thermal inertia):
```math
T_{ss} = \\left[\\frac{(1 - R_{vis}) \\Phi_\\odot}{\\varepsilon \\sigma}\\right]^{1/4}
```
where:
- ``\\Phi_\\odot = \\Phi_0 / r^2`` is the solar flux at distance r [au]
- ``\\Phi_0 = 1366`` W/m² is the solar constant at 1 au
- ``\\sigma`` is the Stefan-Boltzmann constant

# Notes
- This gives the maximum temperature for a non-rotating asteroid
- Actual subsolar temperature may be lower due to thermal inertia
- Valid for airless bodies with negligible heat conduction
"""
subsolar_temperature(r☉, params::AbstractThermoParams) = subsolar_temperature(r☉, params.reflectance_vis, params.emissivity)

function subsolar_temperature(r☉, R_vis, ε)
    Φ = SOLAR_CONST / (norm(r☉) * m2au)^2  # Energy flux at the solar distance [W/m²]
    Tₛₛ = ((1 - R_vis) * Φ / (ε * σ_SB))^(1/4)

    return Tₛₛ
end


"""
    init_temperature!(stpm::SingleAsteroidThermoPhysicalModel, T₀::Real)

Initialize all temperature cells at the given temperature `T₀`

# Arguments
- `stpm` : Thermophysical model for a single asteroid
- `T₀`   : Initial temperature of all cells [K]
"""
function init_temperature!(stpm::SingleAsteroidThermoPhysicalModel, T₀::Real)
    stpm.temperature .= T₀
end


"""
    init_temperature!(btpm::BinaryAsteroidThermoPhysicalModel, T₀::Real)

Initialize all temperature cells at the given temperature `T₀`

# Arguments
- `btpm` : Thermophysical model for a binary asteroid
- `T₀`   : Initial temperature of all cells [K]
"""
function init_temperature!(btpm::BinaryAsteroidThermoPhysicalModel, T₀::Real)
    init_temperature!(btpm.pri, T₀)
    init_temperature!(btpm.sec, T₀)
end


# ╔═══════════════════════════════════════════════════════════════════╗
# ║                    Thermophysical modeling                        ║
# ╚═══════════════════════════════════════════════════════════════════╝


"""
    run_TPM!(stpm::SingleAsteroidThermoPhysicalModel, ephem, times_to_save, face_ID; show_progress=true) -> result

Execute the thermophysical model simulation for a single asteroid over the specified time period.

# Arguments
- `stpm::SingleAsteroidThermoPhysicalModel` : Thermophysical model containing shape, thermal parameters, and state
- `ephem` : Ephemerides data structure containing:
    - `time::Vector{Float64}`   : Time points for the simulation [s]
    - `sun::Vector{SVector{3}}` : Sun's position vectors in the asteroid-fixed frame (not normalized) [m]
- `times_to_save::Vector{Float64}` : Specific time points at which to save detailed temperature data [s]
- `face_ID::Vector{Int}` : Face indices for which to save subsurface temperature profiles

# Keyword Arguments
- `show_progress::Bool=true` : Display progress meter during simulation

# Returns
- `result::SingleAsteroidThermoPhysicalModelResult` : Structure containing:
    - Time series of energy balance (E_in, E_out)
    - Thermal forces and torques at each time step
    - Surface temperatures at specified save times
    - Subsurface temperature profiles for selected faces

# Algorithm
1. For each time step:
   - Update solar flux based on sun position
   - Calculate scattered light flux (if self-heating enabled)
   - Calculate thermal radiation flux (if self-heating enabled)
   - Compute thermal forces and torques
   - Save results if at a save point
   - Update temperature distribution for next step

# Example
```julia
# Setup ephemerides
ephem = (
    time = collect(0:60:86400),  # One day, 1-minute steps
    sun = [SVector(au2m, 0.0, 0.0) for _ in 1:1441]  # Sun at 1 AU
)

# Run simulation
times_to_save = [0.0, 21600.0, 43200.0, 64800.0, 86400.0]  # Every 6 hours
face_ID = [1, 100, 500]  # Save subsurface data for these faces
result = run_TPM!(stpm, ephem, times_to_save, face_ID)

# Check energy balance
println("Final E_out/E_in ratio: ", result.E_out[end]/result.E_in[end])
```

# Performance Notes
- Computation time scales with number of faces and time steps
- Self-shadowing and self-heating calculations add significant overhead

# See Also
- `init_temperature!` to set initial conditions
- `export_TPM_results` to save results to files
- `update_temperature!` for the core temperature update algorithm
"""
function run_TPM!(stpm::SingleAsteroidThermoPhysicalModel, ephem, times_to_save::Vector{Float64}, face_ID::Vector{Int}; show_progress=true)

    result = SingleAsteroidThermoPhysicalModelResult(stpm, ephem, times_to_save, face_ID)

    ## ProgressMeter setting
    if show_progress
        p = Progress(length(ephem.time); dt=1, desc="Running TPM...", showspeed=true)
        ProgressMeter.ijulia_behavior(:clear)
    end

    for i_time in eachindex(ephem.time)
        r☉ = ephem.sun[i_time]

        update_flux_all!(stpm, r☉)                # Update all energy fluxes to the surface
        update_thermal_force!(stpm)               # Calculate thermal forces and torques
        update_TPM_result!(result, stpm, i_time)  # Save data

        ## Update the progress meter
        if show_progress
            showvalues = [
                ("Timestep     ", i_time),
                ("E_out / E_in ", result.E_out[i_time] / result.E_in[i_time]),  # Energy output/input ratio at the time step
            ]
            ProgressMeter.next!(p; showvalues)
        end

        i_time == length(ephem.time) && break  # Stop to update the temperature at the final step
        Δt = ephem.time[i_time+1] - ephem.time[i_time]
        update_temperature!(stpm, Δt)
    end

    return result
end

"""
    run_TPM!(btpm::BinaryAsteroidThermoPhysicalModel, ephem, times_to_save, face_ID_pri, face_ID_sec; show_progress=true) -> result

Execute the thermophysical model simulation for a binary asteroid system over the specified time period.

# Arguments
- `btpm::BinaryAsteroidThermoPhysicalModel` : Thermophysical model for a binary asteroid
- `ephem` : Ephemerides data structure containing:
    - `time::Vector{Float64}`     : Time points for the simulation [s]
    - `sun::Vector{SVector{3}}`   : Sun's position vectors in the primary's frame (NOT normalized) [m]
    - `sec::Vector{SVector{3}}`   : Secondary's position vectors in the primary's frame [m]
    - `P2S::Vector{SMatrix{3,3}}` : Rotation matrices from primary to secondary frames
- `times_to_save::Vector{Float64}` : Specific time points at which to save detailed temperature data [s]
- `face_ID_pri::Vector{Int}` : Face indices for which to save subsurface temperature profiles for the primary
- `face_ID_sec::Vector{Int}` : Face indices for which to save subsurface temperature profiles for the secondary

# Keyword Arguments
- `show_progress::Bool=true` : Display progress meter during simulation

# Returns
- `result::BinaryAsteroidThermoPhysicalModelResult` : Structure containing results for both components

# Algorithm
1. For each time step:
   - Update all energy fluxes using the unified API (solar+eclipse, scattered, thermal radiation, mutual heating)
   - Compute thermal forces and torques for both components
   - Save results if at a save point
   - Update temperature distribution for next time step

# Notes
- Uses the unified `update_flux_all!` API which internally handles all coordinate transformations
- Automatically respects all shadowing and heating flags (SELF_SHADOWING, SELF_HEATING, MUTUAL_SHADOWING, MUTUAL_HEATING)
- Compatible with AsteroidShapeModels.jl v0.4.1 which fixes critical eclipse shadowing bugs

# See Also
- `init_temperature!` to set initial conditions
- `export_TPM_results` to save results to files
- `update_flux_all!` for the unified flux update API
"""
function run_TPM!(btpm::BinaryAsteroidThermoPhysicalModel, ephem, times_to_save::Vector{Float64}, face_ID_pri::Vector{Int}, face_ID_sec::Vector{Int}; show_progress=true)

    result = BinaryAsteroidThermoPhysicalModelResult(btpm, ephem, times_to_save, face_ID_pri, face_ID_sec)

    ## ProgressMeter setting
    if show_progress
        p = Progress(length(ephem.time); dt=1, desc="Running TPM...", showspeed=true)
        ProgressMeter.ijulia_behavior(:clear)
    end

    for i_time in eachindex(ephem.time)
        r☉₁ = ephem.sun[i_time]  # Sun's position in the primary's frame
        r₁₂ = ephem.sec[i_time]  # Secondary's position in the primary's frame
        R₁₂ = ephem.P2S[i_time]  # Rotation matrix from primary to secondary frames

        update_flux_all!(btpm, r☉₁, r₁₂, R₁₂)     # Update all energy fluxes to surface
        update_thermal_force!(btpm)               # Calculate thermal forces and torques
        update_TPM_result!(result, btpm, i_time)  # Save data

        ## Update the progress meter
        if show_progress
            showvalues = [
                ("Timestep                   ", i_time),
                ("E_out / E_in for primary   ", result.pri.E_out[i_time] / result.pri.E_in[i_time]),  # Energy output/input ratio on the primary at the time step
                ("E_out / E_in for secondary ", result.sec.E_out[i_time] / result.sec.E_in[i_time]),  # Energy output/input ratio on the secondary at the time step
            ]
            ProgressMeter.next!(p; showvalues)
        end

        ## Update temperature distribution
        i_time == length(ephem.time) && break  # Stop to update the temperature at the final step
        Δt = ephem.time[i_time+1] - ephem.time[i_time]
        update_temperature!(btpm, Δt)
    end

    return result
end
