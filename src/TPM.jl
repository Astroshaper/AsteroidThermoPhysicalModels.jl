
# ****************************************************************
#         Types of solvers for heat conduction equations
# ****************************************************************

"""
Abstract type of a solver for a heat conduction equation 
"""
abstract type HeatConductionSolver end


"""
Type of the forward Euler method:
- Explicit in time
- First order in time

The `ForwardEulerSolver` type includes a vector for the temperature at the next time step.
"""
struct ForwardEulerSolver <: HeatConductionSolver
    T::Vector{Float64}
end

ForwardEulerSolver(thermo_params::AbstractThermoParams) = ForwardEulerSolver(thermo_params.Nz)
ForwardEulerSolver(N::Integer) = ForwardEulerSolver(zeros(N))


"""
Type of the backward Euler method:
- Implicit in time (Unconditionally stable in the heat conduction equation)
- First order in time

The `BackwardEulerSolver` type has vectors for the tridiagonal matrix algorithm.
"""
struct BackwardEulerSolver <: HeatConductionSolver
    a::Vector{Float64}
    b::Vector{Float64}
    c::Vector{Float64}
    d::Vector{Float64}
    x::Vector{Float64}
end

BackwardEulerSolver(thermo_params::AbstractThermoParams) = BackwardEulerSolver(thermo_params.Nz)
BackwardEulerSolver(N::Integer) = BackwardEulerSolver(zeros(N), zeros(N), zeros(N), zeros(N), zeros(N))


"""
Type of the Crank-Nicolson method:
- Implicit in time (Unconditionally stable in the heat conduction equation)
- Second order in time

The `CrankNicolsonSolver` type has vectors for the tridiagonal matrix algorithm.

# References
- https://en.wikipedia.org/wiki/Crank–Nicolson_method
"""
struct CrankNicolsonSolver <: HeatConductionSolver
    a::Vector{Float64}
    b::Vector{Float64}
    c::Vector{Float64}
    d::Vector{Float64}
    x::Vector{Float64}
end

CrankNicolsonSolver(thermo_params::AbstractThermoParams) = CrankNicolsonSolver(thermo_params.Nz)
CrankNicolsonSolver(N::Integer) = CrankNicolsonSolver(zeros(N), zeros(N), zeros(N), zeros(N), zeros(N))


# ****************************************************************
#    Types of boundary conditions for heat conduction equations
# ****************************************************************

"""
Abstract type of a boundary condition for a heat conduction equation
"""
abstract type BoundaryCondition end


"""
Singleton type of radiation boundary condition
"""
struct RadiationBoundaryCondition <: BoundaryCondition end


"""
Singleton type of insulation boundary condition
"""
struct InsulationBoundaryCondition <: BoundaryCondition end


"""
Type of isothermal boundary condition
"""
struct IsothermalBoundaryCondition <: BoundaryCondition
    T_iso::Float64
end


# ****************************************************************
#                Types of thermophsycal models
# ****************************************************************


"""
Abstract type of a thermophysical model
"""
abstract type ThermoPhysicalModel end


"""
    struct SingleTPM <: ThermoPhysicalModel

# Fields
- `shape`          : Shape model
- `thermo_params`  : Thermophysical parameters

- `flux`           : Flux on each face. Matrix of size (Number of faces, 3). Three components are:
    - `flux[:, 1]`     : F_sun,  flux of direct sunlight
    - `flux[:, 2]`     : F_scat, flux of scattered light
    - `flux[:, 3]`     : F_rad,  flux of thermal emission from surrounding surface
- `temperature`    : Temperature matrix `(Nz, Ns)` according to the number of depth cells `Nz` and the number of faces `Ns`.

- `face_forces`    : Thermal force on each face
- `force`          : Thermal recoil force at body-fixed frame (Yarkovsky effect)
- `torque`         : Thermal recoil torque at body-fixed frame (YORP effect)

- `SELF_SHADOWING` : Flag to consider self-shadowing
- `SELF_HEATING`   : Flag to consider self-heating
- `SOLVER`         : Solver of heat conduction equation
- `BC_UPPER`       : Boundary condition at the upper boundary
- `BC_LOWER`       : Boundary condition at the lower boundary

# TO DO:
- roughness_maps   ::ShapeModel[]
"""
struct SingleTPM <: ThermoPhysicalModel
    shape          ::ShapeModel
    thermo_params  ::Union{UniformThermoParams, NonUniformThermoParams}

    flux           ::Matrix{Float64}  # (Ns, 3)
    temperature    ::Matrix{Float64}  # (Nz, Ns)

    face_forces    ::Vector{SVector{3, Float64}}
    force          ::MVector{3, Float64}
    torque         ::MVector{3, Float64}

    SELF_SHADOWING ::Bool
    SELF_HEATING   ::Bool
    SOLVER         ::Union{ForwardEulerSolver, BackwardEulerSolver, CrankNicolsonSolver}
    BC_UPPER       ::Union{RadiationBoundaryCondition, InsulationBoundaryCondition, IsothermalBoundaryCondition}
    BC_LOWER       ::Union{RadiationBoundaryCondition, InsulationBoundaryCondition, IsothermalBoundaryCondition}
end


"""
    SingleTPM(shape, thermo_params; SELF_SHADOWING=true, SELF_HEATING=true) -> stpm

Construct a thermophysical model for a single asteroid (`SingleTPM`).

# Arguments
- `shape`          : Shape model
- `thermo_params`  : Thermophysical parameters

# Keyword arguments
- `SELF_SHADOWING` : Flag to consider self-shadowing
- `SELF_HEATING`   : Flag to consider self-heating
- `SOLVER`         : Solver of heat conduction equation
- `BC_UPPER`       : Boundary condition at the upper boundary
- `BC_LOWER`       : Boundary condition at the lower boundary
"""
function SingleTPM(shape, thermo_params; SELF_SHADOWING, SELF_HEATING, SOLVER, BC_UPPER, BC_LOWER)

    Nz = thermo_params.Nz
    Ns = length(shape.faces)

    flux = zeros(Ns, 3)
    temperature = zeros(Nz, Ns)

    face_forces = [zero(SVector{3, Float64}) for _ in shape.faces]
    force  = zero(MVector{3, Float64})
    torque = zero(MVector{3, Float64})

    SingleTPM(shape, thermo_params, flux, temperature, face_forces, force, torque, SELF_SHADOWING, SELF_HEATING, SOLVER, BC_UPPER, BC_LOWER)
end


"""
    struct BinaryTPM <: ThermoPhysicalModel

# Fields
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


"""
    BinaryTPM(pri, sec; MUTUAL_SHADOWING=true, MUTUAL_HEATING=true) -> btpm

Construct a thermophysical model for a binary asteroid (`BinaryTPM`).
"""
function BinaryTPM(pri, sec; MUTUAL_SHADOWING, MUTUAL_HEATING)
    BinaryTPM(pri, sec, MUTUAL_SHADOWING, MUTUAL_HEATING)
end


"""
Return surface temperature of a single asteroid corrsponding to each face.
"""
surface_temperature(stpm::SingleTPM) = stpm.temperature[begin, :]


# ****************************************************************
#                     Output data format
# ****************************************************************

"""
    struct SingleTPMResult

Output data format for `SingleTPM` 

# Fields
## Saved at all time steps
- `E_in`   : Input energy per second on the whole surface [W]
- `E_out`  : Output enegey per second from the whole surface [W]
- `force`  : Thermal force on the asteroid [N]
- `torque` : Thermal torque on the asteroid [N ⋅ m]

## Saved only at the time steps desired by the user
- `time_begin` : Time to start storing temperature
- `time_end`   : Time to finish storing temperature
- `face_id`    : Indices of faces at which you want to save temperature distribution in depth direction
- `surf_temp`  : Surface temperature, a matrix in size of `(Ns, Nt)`.
    - `Ns` : Number of faces
    - `Nt` : Number of time steps to save surface temperature
- `face_temp`  : Temperature as a function of depth and time, an array in size of `(Nz, Ns, Nt)`.
    - `Nz` : The number of the depth nodes
    - `Ns` : The number of faces to save temperature
    - `Nt` : The number of time steps to save temperature
"""
struct SingleTPMResult
    E_in       ::Vector{Float64}
    E_out      ::Vector{Float64}
    force      ::Vector{SVector{3, Float64}}
    torque     ::Vector{SVector{3, Float64}}

    time_begin ::Float64
    time_end   ::Float64
    face_id    ::Vector{Int}
    surf_temp  ::Matrix{Float64}
    face_temp  ::Array{Float64, 3}
end


"""
Outer constructor of `SingleTPMResult`

# Arguments
- `stpm`    : Thermophysical model for a single asteroid
- `ephem`   : Ephemerides
- `time_id` : Indices of time steps at which you want to save temperature
- `face_id` : Indices of faces at which you want to save temperature distribution in depth direction
"""
function SingleTPMResult(stpm::SingleTPM, ephem, time_begin::Real, time_end::Real, face_id::Vector{Int})
    E_in   = zeros(length(ephem.time))
    E_out  = zeros(length(ephem.time))
    force  = zeros(SVector{3, Float64}, length(ephem.time))
    torque = zeros(SVector{3, Float64}, length(ephem.time))

    Nt_save = count(@. time_begin ≤ ephem.time < time_end)  # Number of time steps to save temperature

    surf_temp = zeros(length(stpm.shape.faces), Nt_save)
    face_temp = zeros(stpm.thermo_params.Nz, length(face_id), Nt_save)    

    return SingleTPMResult(E_in, E_out, force, torque, time_begin, time_end, face_id, surf_temp, face_temp)
end


"""
    struct BinaryTPMResult

Output data format for `BinaryTPM`

# Fields
- `pri` : TPM result for the primary
- `sec` : TPM result for the secondary
"""
struct BinaryTPMResult
    pri::SingleTPMResult
    sec::SingleTPMResult
end


"""
Outer constructor of `BinaryTPMResult`

# Arguments
- `btpm`    : Thermophysical model for a binary asteroid
"""
function BinaryTPMResult(btpm::BinaryTPM)
    # pri = SingleTPMResult(btpm.pri, btpm.pri.ephem, btpm.pri.time_begin, btpm.pri.time_end, btpm.pri.face_id)
    # sec = SingleTPMResult(btpm.sec, btpm.sec.ephem, btpm.sec.time_begin, btpm.sec.time_end, btpm.sec.face_id)

    # return BinaryTPMResult(pri, sec)
end



function save_TPM_result!(result::SingleTPMResult, stpm::SingleTPM, ephem, nₜ::Integer)
    result.E_in[nₜ]   = energy_in(stpm)
    result.E_out[nₜ]  = energy_out(stpm)
    result.force[nₜ]  = stpm.force
    result.torque[nₜ] = stpm.torque

    if result.time_begin ≤ ephem.time[nₜ] < result.time_end   # if you want to save temperature at this time step
        nₜ_offset = count(@. ephem.time < result.time_begin)  # Index-offset before storing temperature
        nₜ_save = nₜ - nₜ_offset
        result.surf_temp[:, nₜ_save] .= surface_temperature(stpm)

        for (nₛ_save, nₛ) in enumerate(result.face_id)
            result.face_temp[:, nₛ_save, nₜ_save] .= stpm.temperature[:, nₛ]
        end
    end
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
    init_temperature!(stpm::SingleTPM, T₀::Real)

Initialize all temperature cells at the given temperature `T₀`

# Arguments
- `stpm` : Thermophysical model for a single asteroid
- `T₀`   : Initial temperature of all cells [K]
"""
function init_temperature!(stpm::SingleTPM, T₀::Real)
    stpm.temperature .= T₀
end


"""
    init_temperature!(btpm::BinaryTPM, T₀::Real)

Initialize all temperature cells at the given temperature `T₀`

# Arguments
- `btpm` : Thermophysical model for a binary asteroid
- `T₀`   : Initial temperature of all cells [K]
"""
function init_temperature!(btpm::BinaryTPM, T₀::Real)
    init_temperature!(btpm.pri, T₀)
    init_temperature!(btpm.sec, T₀)
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
function run_TPM!(stpm::SingleTPM, ephem, time_begin::Real, time_end::Real, face_id::Vector{Int})

    result = SingleTPMResult(stpm, ephem, time_begin, time_end, face_id)

    ## ProgressMeter setting
    p = Progress(length(ephem.time); dt=1, desc="Running TPM...", showspeed=true)
    ProgressMeter.ijulia_behavior(:clear)
    
    @time for nₜ in eachindex(ephem.time)
        r☉ = ephem.sun[nₜ]

        update_flux_sun!(stpm, r☉)
        update_flux_scat_single!(stpm)
        update_flux_rad_single!(stpm)
        
        update_thermal_force!(stpm)

        save_TPM_result!(result, stpm, ephem, nₜ)  # Save data
        
        ## Update the progress meter
        showvalues = [("Timestep", nₜ), ("E_out / E_in", result.E_out[nₜ] / result.E_in[nₜ])]
        ProgressMeter.next!(p; showvalues)

        nₜ == length(ephem.time) && break  # Stop to update the temperature at the final step
        Δt = ephem.time[nₜ+1] - ephem.time[nₜ]
        update_temperature!(stpm, Δt)
    end

    return result
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

    @time for nₜ in eachindex(ephem.time)
        r☉₁ = ephem.sun1[nₜ]  # Sun's position in the primary's frame
        r☉₂ = ephem.sun2[nₜ]  # Sun's position in the secondary's frame
        rₛ  = ephem.sec[nₜ]   # Secondary's position in the primary's frame
        R₂₁ = ephem.S2P[nₜ]   # Rotation matrix from secondary to primary frames

        ## Update enegey flux
        update_flux_sun!(btpm, r☉₁, r☉₂)
        update_flux_scat_single!(btpm)
        update_flux_rad_single!(btpm)

        mutual_shadowing!(btpm, r☉₁, rₛ, R₂₁)  # Mutual-shadowing (eclipse)
        mutual_heating!(btpm, rₛ, R₂₁)         # Mutual-heating

        update_thermal_force!(btpm)

        ## Save data for primary
        surf_temps[1][:, nₜ] .= surface_temperature(btpm.pri)
        forces[1][nₜ]  .= btpm.pri.force   # Body-fixed frame
        torques[1][nₜ] .= btpm.pri.torque  # Body-fixed frame

        ## Save data for secondary
        surf_temps[2][:, nₜ] .= surface_temperature(btpm.sec)
        forces[2][nₜ]  .= btpm.sec.force   # Body-fixed frame
        torques[2][nₜ] .= btpm.sec.torque  # Body-fixed frame
    
        ## Energy input/output
        E_cons_pri = energy_out(btpm.pri) / energy_in(btpm.pri)
        E_cons_sec = energy_out(btpm.sec) / energy_in(btpm.sec)

        ## Update the progress meter
        showvalues = [("Timestep", nₜ), ("E_cons for primary", E_cons_pri), ("E_cons for secondary", E_cons_sec)]
        ProgressMeter.next!(p; showvalues)
        
        ## Update temperature distribution
        nₜ == length(ephem.time) && break  # Stop to update the temperature at the final step
        Δt = ephem.time[nₜ+1] - ephem.time[nₜ]
        update_temperature!(btpm, Δt)
    end
    
    jldsave(savepath; btpm, ephem, surf_temps, forces, torques)
end

