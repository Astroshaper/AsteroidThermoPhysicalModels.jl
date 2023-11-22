
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
- `times`  : Timesteps, given the same vector as `ephem.time` [s]
- `E_in`   : Input energy per second on the whole surface [W]
- `E_out`  : Output enegey per second from the whole surface [W]
- `E_cons` : Energy conservation ratio [-], ratio of total energy going out to total energy coming in in the last rotation cycle
- `force`  : Thermal force on the asteroid [N]
- `torque` : Thermal torque on the asteroid [N ⋅ m]

## Saved only at the time steps desired by the user
- `times_to_save` : Timesteps to save temperature [s]
- `surf_temp`     : Surface temperature [K], a matrix in size of `(Ns, Nt)`.
    - `Ns` : Number of faces
    - `Nt` : Number of time steps to save surface temperature
- `face_temp`     : Temperature [K] as a function of depth [m] and time [s], `Dict` with face ID as key and a matrix `(Nz, Nt)` as an entry.
    - `Nz` : The number of the depth nodes
    - `Nt` : The number of time steps to save temperature
"""
struct SingleTPMResult
    times  ::Vector{Float64}
    E_in   ::Vector{Float64}
    E_out  ::Vector{Float64}
    E_cons ::Vector{Union{Float64, Missing}}
    force  ::Vector{SVector{3, Float64}}
    torque ::Vector{SVector{3, Float64}}

    times_to_save ::Vector{Float64}
    surf_temp     ::Matrix{Float64}
    face_temp     ::Dict{Int, Matrix{Float64}}
end


"""
Outer constructor of `SingleTPMResult`

# Arguments
- `stpm`          : Thermophysical model for a single asteroid
- `ephem`         : Ephemerides
- `times_to_save` : Timesteps to save temperature
- `face_ID`       : Face indices where to save subsurface temperature
"""
function SingleTPMResult(stpm::SingleTPM, ephem, times_to_save::Vector{Float64}, face_ID::Vector{Int})
    E_in   = zeros(length(ephem.time))
    E_out  = zeros(length(ephem.time))
    E_cons = Vector{Union{Float64, Missing}}(missing, length(ephem.time))
    force  = zeros(SVector{3, Float64}, length(ephem.time))
    torque = zeros(SVector{3, Float64}, length(ephem.time))

    surf_temp = zeros(length(stpm.shape.faces), length(times_to_save))
    face_temp = Dict{Int, Matrix{Float64}}(
        nₛ => zeros(stpm.thermo_params.Nz, length(times_to_save)) for nₛ in face_ID
    )

    return SingleTPMResult(ephem.time, E_in, E_out, E_cons, force, torque, times_to_save, surf_temp, face_temp)
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
- `btpm`          : Thermophysical model for a binary asteroid
- `ephem`         : Ephemerides
- `times_to_save` : Timesteps to save temperature (Common to both the primary and the secondary)
- `face_ID_pri`   : Face indices at which you want to save temperature as a function of depth for the primary
- `face_ID_sec`   : Face indices at which you want to save temperature as a function of depth for the secondary
"""
function BinaryTPMResult(btpm::BinaryTPM, ephem, times_to_save::Vector{Float64}, face_ID_pri::Vector{Int}, face_ID_sec::Vector{Int})
    result_pri = SingleTPMResult(btpm.pri, ephem, times_to_save, face_ID_pri)
    result_sec = SingleTPMResult(btpm.sec, ephem, times_to_save, face_ID_sec)

    return BinaryTPMResult(result_pri, result_sec)
end


"""
    update_TPM_result!(result::SingleTPMResult, stpm::SingleTPM, ephem, nₜ::Integer)

Save the results of TPM at the time step `nₜ` to `result`.

# Arguments
- `result` : Output data format for `SingleTPM`
- `stpm`   : Thermophysical model for a single asteroid
- `ephem`  : Ephemerides
- `nₜ`     : Time step to save data
"""
function update_TPM_result!(result::SingleTPMResult, stpm::SingleTPM, nₜ::Integer)
    result.E_in[nₜ]   = energy_in(stpm)
    result.E_out[nₜ]  = energy_out(stpm)
    result.force[nₜ]  = stpm.force
    result.torque[nₜ] = stpm.torque

    P  = stpm.thermo_params.P  # Rotation period
    t  = result.times[nₜ]        # Current time
    t₀ = result.times[begin]     # Time at the beginning of the simulation

    if t > t₀ + P  # Note that `E_cons` cannot be calculated during the first rotation
        nsteps_in_period = count(@. t - P ≤ result.times < t)  # Number of time steps within the last rotation

        ΣE_in  = sum(result.E_in[n-1]  * (result.times[n] - result.times[n-1]) for n in (nₜ - nsteps_in_period + 1):nₜ)
        ΣE_out = sum(result.E_out[n-1] * (result.times[n] - result.times[n-1]) for n in (nₜ - nsteps_in_period + 1):nₜ)

        result.E_cons[nₜ] = ΣE_out / ΣE_in
    end

    if t in result.times_to_save  # In the step of saving temperature
        nₜ_save = findfirst(isequal(t), result.times_to_save)

        result.surf_temp[:, nₜ_save] .= surface_temperature(stpm)

        for (nₛ, temp) in result.face_temp
            temp[:, nₜ_save] .= stpm.temperature[:, nₛ]
        end
    end
end


"""
    update_TPM_result!(result::BinaryTPMResult, btpm::BinaryTPM, ephem, nₜ::Integer)

Save the results of TPM at the time step `nₜ` to `result`.

# Arguments
- `result` : Output data format for `BinaryTPM`
- `btpm`   : Thermophysical model for a binary asteroid
- `ephem`  : Ephemerides
- `nₜ`     : Time step
"""
function update_TPM_result!(result::BinaryTPMResult, btpm::BinaryTPM, nₜ::Integer)
    update_TPM_result!(result.pri, btpm.pri, nₜ)
    update_TPM_result!(result.sec, btpm.sec, nₜ)
end


"""
    export_TPM_results(dirpath, result::SingleTPMResult)

Export the result of `SingleTPM` to CSV files.

# Arguments
- `dirpath` :  Path to the directory to save CSV files
- `result`  : Output data format for `SingleTPM`

# TO DO
- Save the depths of the calculation nodes
- Save README for the data file
"""
function export_TPM_results(dirpath, result::SingleTPMResult)
    
    df = DataFrame()
    df.time     = result.times
    df.E_in     = result.E_in
    df.E_out    = result.E_out
    df.E_cons   = result.E_cons
    df.force_x  = [f[1] for f in result.force]
    df.force_y  = [f[2] for f in result.force]
    df.force_z  = [f[3] for f in result.force]
    df.torque_x = [τ[1] for τ in result.torque]
    df.torque_y = [τ[2] for τ in result.torque]
    df.torque_z = [τ[3] for τ in result.torque]
    
    CSV.write(joinpath(dirpath, "data.csv"), df)

    ##= Surface temperature =##
    header = string.(result.times_to_save)
    CSV.write(
        joinpath(dirpath, "surf_temp.csv"),
        DataFrame(result.surf_temp, header)
    )

    ##= Subsurface temperature =##
    for (nₛ, temp) in result.face_temp
        CSV.write(
            joinpath(dirpath, "face_temp_$(lpad(nₛ, 7, '0')).csv"),
            DataFrame(temp, header)
        )
    end
end


"""
    export_TPM_results(filepath, result::BinaryTPMResult)

Export the result of `BinaryTPM` to CSV files.

# Arguments
- `dirpath` : Path to the directory to save CSV files
- `result`  : Output data format for `BinaryTPM`
"""
function export_TPM_results(dirpath, result::BinaryTPMResult)
    dirpath_pri = joinpath(dirpath, "pri")
    dirpath_sec = joinpath(dirpath, "sec")

    mkpath(dirpath_pri)
    mkpath(dirpath_sec)

    export_TPM_results(dirpath_pri, result.pri)
    export_TPM_results(dirpath_sec, result.sec)
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
- `stpm`          : Thermophysical model for a single asteroid
- `ephem`         : Ephemerides
    - `ephem.time` : Ephemeris times
    - `ephem.sun`  : Sun's position in the asteroid-fixed frame (Not normalized)
- `times_to_save` : Timesteps to save temperature
- `face_ID`       : Face indices where to save subsurface termperature

# Keyword arguments
- `show_progress` : Flag to show the progress meter
"""
function run_TPM!(stpm::SingleTPM, ephem, times_to_save::Vector{Float64}, face_ID::Vector{Int}; show_progress=true)

    result = SingleTPMResult(stpm, ephem, times_to_save, face_ID)

    ## ProgressMeter setting
    if show_progress
        p = Progress(length(ephem.time); dt=1, desc="Running TPM...", showspeed=true)
        ProgressMeter.ijulia_behavior(:clear)
    end
    
    for nₜ in eachindex(ephem.time)
        r☉ = ephem.sun[nₜ]

        update_flux_sun!(stpm, r☉)
        update_flux_scat_single!(stpm)
        update_flux_rad_single!(stpm)
        
        update_thermal_force!(stpm)

        update_TPM_result!(result, stpm, nₜ)  # Save data
        
        ## Update the progress meter
        if show_progress
            showvalues = [
                ("Timestep ", nₜ),
                ("E_cons   ", result.E_cons[nₜ]),
            ]
            ProgressMeter.next!(p; showvalues)
        end

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
- `btpm`          : Thermophysical model for a binary asteroid
- `ephem`         : Ephemerides
    - `time` : Ephemeris times
    - `sun1` : Sun's position in the primary's frame
    - `sun2` : Sun's position in the secondary's frame
    - `sec`  : Secondary's position in the primary's frame
    - `P2S`  : Rotation matrix from primary to secondary frames
    - `S2P`  : Rotation matrix from secondary to primary frames
- `times_to_save` : Timesteps to save temperature
- `face_ID_pri`   : Face indices where to save subsurface termperature for the primary
- `face_ID_sec`   : Face indices where to save subsurface termperature for the secondary

# Keyword arguments
- `show_progress` : Flag to show the progress meter
"""
function run_TPM!(btpm::BinaryTPM, ephem, times_to_save::Vector{Float64}, face_ID_pri::Vector{Int}, face_ID_sec::Vector{Int}; show_progress=true)

    result = BinaryTPMResult(btpm, ephem, times_to_save, face_ID_pri, face_ID_sec)

    ## ProgressMeter setting
    if show_progress
        p = Progress(length(ephem.time); dt=1, desc="Running TPM...", showspeed=true)
        ProgressMeter.ijulia_behavior(:clear)
    end
    
    for nₜ in eachindex(ephem.time)
        r☉₁ = ephem.sun1[nₜ]  # Sun's position in the primary's frame
        r☉₂ = ephem.sun2[nₜ]  # Sun's position in the secondary's frame
        rₛ  = ephem.sec[nₜ]   # Secondary's position in the primary's frame
        R₂₁ = ephem.S2P[nₜ]   # Rotation matrix from secondary to primary frames

        ## Update enegey flux
        update_flux_sun!(btpm, r☉₁, r☉₂)
        mutual_shadowing!(btpm, r☉₁, rₛ, R₂₁)  # Mutual-shadowing (eclipse)
        update_flux_scat_single!(btpm)
        update_flux_rad_single!(btpm)
        mutual_heating!(btpm, rₛ, R₂₁)         # Mutual-heating

        update_thermal_force!(btpm)

        update_TPM_result!(result, btpm, nₜ)  # Save data

        ## Update the progress meter
        if show_progress
            showvalues = [
                ("Timestep             ", nₜ),
                ("E_cons for primary   ", result.pri.E_cons[nₜ]),
                ("E_cons for secondary ", result.sec.E_cons[nₜ]),
            ]
            ProgressMeter.next!(p; showvalues)
        end
        
        ## Update temperature distribution
        nₜ == length(ephem.time) && break  # Stop to update the temperature at the final step
        Δt = ephem.time[nₜ+1] - ephem.time[nₜ]
        update_temperature!(btpm, Δt)
    end

    return result
end

