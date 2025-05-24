
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

ForwardEulerSolver(thermo_params::AbstractThermoParams) = ForwardEulerSolver(thermo_params.n_depth)
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

BackwardEulerSolver(thermo_params::AbstractThermoParams) = BackwardEulerSolver(thermo_params.n_depth)
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

CrankNicolsonSolver(thermo_params::AbstractThermoParams) = CrankNicolsonSolver(thermo_params.n_depth)
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
Abstract type of an asteroid's thermophysical model.
The `AbstractAsteroidTPM` type is an alias for `AbstractAsteroidThermoPhysicalModel`.
"""
abstract type AbstractAsteroidThermoPhysicalModel end


"""
    broadcast_thermo_params!(thermo_params::ThermoParams, n_face::Int)

Broadcast the thermophysical parameters to all faces if the values are uniform globally.

# Arguments
- `thermo_params` : Thermophysical parameters
- `n_face`        : Number of faces on the shape model
"""
function broadcast_thermo_params!(thermo_params::ThermoParams, n_face::Int)
    if length(thermo_params.thermal_conductivity) == 1
        resize!(thermo_params.thermal_conductivity, n_face)
        resize!(thermo_params.density,              n_face)
        resize!(thermo_params.heat_capacity,        n_face)
        resize!(thermo_params.reflectance_vis,      n_face)
        resize!(thermo_params.reflectance_ir,       n_face)
        resize!(thermo_params.emissivity,           n_face)

        thermo_params.thermal_conductivity[2:end] .= thermo_params.thermal_conductivity[begin]
        thermo_params.density[2:end]              .= thermo_params.density[begin]
        thermo_params.heat_capacity[2:end]        .= thermo_params.heat_capacity[begin]
        thermo_params.reflectance_vis[2:end]      .= thermo_params.reflectance_vis[begin]
        thermo_params.reflectance_ir[2:end]       .= thermo_params.reflectance_ir[begin]
        thermo_params.emissivity[2:end]           .= thermo_params.emissivity[begin]
    elseif length(thermo_params.thermal_conductivity) == n_face
        # Do nothing
    else
        throw(ArgumentError("The length of the thermophysical parameters is invalid."))
    end
end


"""
    broadcast_thermo_params!(thermo_params::ThermoParams, shape::ShapeModel)

Broadcast the thermophysical parameters to all faces if the values are uniform globally.

# Arguments
- `thermo_params` : Thermophysical parameters
- `shape`         : Shape model
"""
broadcast_thermo_params!(thermo_params::ThermoParams, shape::ShapeModel) = broadcast_thermo_params!(thermo_params, length(shape.faces))


"""
    struct SingleAsteroidThermoPhysicalModel <: AbstractAsteroidThermoPhysicalModel

# Fields
- `shape`          : Shape model
- `thermo_params`  : Thermophysical parameters

- `flux_sun`       : Flux of direct sunlight on each face [W/m²]
- `flux_scat`      : Flux of scattered light on each face [W/m²]
- `flux_rad`       : Flux of thermal emission from surrounding surface on each face [W/m²]
- `temperature`    : Temperature matrix `(n_depth, n_face)` according to the number of depth cells `n_depth` and the number of faces `n_face`.

- `face_forces`    : Thermal force on each face
- `force`          : Thermal recoil force at body-fixed frame (Yarkovsky effect)
- `torque`         : Thermal recoil torque at body-fixed frame (YORP effect)

- `SELF_SHADOWING` : Flag to consider self-shadowing
- `SELF_HEATING`   : Flag to consider self-heating
- `SOLVER`         : Solver of heat conduction equation
- `BC_UPPER`       : Boundary condition at the upper boundary
- `BC_LOWER`       : Boundary condition at the lower boundary

# TODO:
- roughness_maps   ::ShapeModel[]
"""
struct SingleAsteroidThermoPhysicalModel{P<:AbstractThermoParams, S<:HeatConductionSolver, BU<:BoundaryCondition, BL<:BoundaryCondition} <: AbstractAsteroidThermoPhysicalModel
    shape          ::ShapeModel
    thermo_params  ::P

    flux_sun       ::Vector{Float64}  # Flux of direct sunlight
    flux_scat      ::Vector{Float64}  # Flux of scattered light
    flux_rad       ::Vector{Float64}  # Flux of thermal emission from surrounding surface
    temperature    ::Matrix{Float64}  # (n_depth, n_face)

    face_forces    ::Vector{SVector{3, Float64}}
    force          ::MVector{3, Float64}
    torque         ::MVector{3, Float64}

    SELF_SHADOWING ::Bool
    SELF_HEATING   ::Bool
    SOLVER         ::S
    BC_UPPER       ::BU
    BC_LOWER       ::BL
end


"""
    SingleAsteroidThermoPhysicalModel(shape, thermo_params; SELF_SHADOWING=true, SELF_HEATING=true) -> stpm

Construct a thermophysical model for a single asteroid (`SingleAsteroidThermoPhysicalModel`).

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
function SingleAsteroidThermoPhysicalModel(shape, thermo_params; SELF_SHADOWING, SELF_HEATING, SOLVER, BC_UPPER, BC_LOWER)

    broadcast_thermo_params!(thermo_params, shape)

    broadcast_thermo_params!(thermo_params, shape)

    n_depth = thermo_params.n_depth
    n_face = length(shape.faces)

    flux_sun = zeros(n_face)
    flux_scat = zeros(n_face)
    flux_rad = zeros(n_face)
    temperature = zeros(n_depth, n_face)

    face_forces = zeros(SVector{3, Float64}, n_face)
    force  = zero(MVector{3, Float64})
    torque = zero(MVector{3, Float64})

    SingleAsteroidThermoPhysicalModel(shape, thermo_params, flux_sun, flux_scat, flux_rad, temperature, face_forces, force, torque, SELF_SHADOWING, SELF_HEATING, SOLVER, BC_UPPER, BC_LOWER)
end


"""
    struct BinaryAsteroidThermoPhysicalModel{M1, M2} <: AbstractAsteroidThermoPhysicalModel

# Fields
- `pri`              : TPM for the primary
- `sec`              : TPM for the secondary
- `MUTUAL_SHADOWING` : Flag to consider mutual shadowing
- `MUTUAL_HEATING`   : Flag to consider mutual heating
"""
struct BinaryAsteroidThermoPhysicalModel{M1, M2} <: AbstractAsteroidThermoPhysicalModel
    pri              ::M1
    sec              ::M2

    MUTUAL_SHADOWING ::Bool
    MUTUAL_HEATING   ::Bool
end


"""
    BinaryAsteroidThermoPhysicalModel(pri, sec; MUTUAL_SHADOWING=true, MUTUAL_HEATING=true) -> btpm

Construct a thermophysical model for a binary asteroid (`BinaryAsteroidThermoPhysicalModel`).
"""
function BinaryAsteroidThermoPhysicalModel(pri, sec; MUTUAL_SHADOWING, MUTUAL_HEATING)
    BinaryAsteroidThermoPhysicalModel(pri, sec, MUTUAL_SHADOWING, MUTUAL_HEATING)
end


"""
Return surface temperature of a single asteroid corrsponding to each face.
"""
surface_temperature(stpm::SingleAsteroidThermoPhysicalModel) = stpm.temperature[begin, :]


# ****************************************************************
#                     Output data format
# ****************************************************************

"""
    struct SingleAsteroidThermoPhysicalModelResult

Output data format for `SingleAsteroidThermoPhysicalModel` 

# Fields
## Saved at all time steps
- `times`  : Timesteps, given the same vector as `ephem.time` [s]
- `E_in`   : Input energy per second on the whole surface [W]
- `E_out`  : Output energey per second from the whole surface [W]
- `force`  : Thermal force on the asteroid [N]
- `torque` : Thermal torque on the asteroid [N ⋅ m]

## Saved only at the time steps desired by the user
- `times_to_save`          : Timesteps to save temperature and thermal force on every face [s]
- `depth_nodes`            : Depths of the calculation nodes for 1-D heat conduction [m], a vector of size `n_depth`
- `surface_temperature`    : Surface temperature [K], a matrix in size of `(n_face, n_time)`.
    - `n_face` : Number of faces
    - `n_time` : Number of time steps to save surface temperature
- `subsurface_temperature` : Temperature [K] as a function of depth [m] and time [s], `Dict` with face ID as key and a matrix `(n_depth, n_time)` as an entry.
    - `n_depth` : The number of the depth nodes
    - `n_time` : The number of time steps to save temperature
- `face_forces`            : Thermal force on every face of the shape model [N], a matrix in size of `(n_face, n_time)`.
    - `n_face` : Number of faces
    - `n_time` : Number of time steps to save surface temperature
"""
struct SingleAsteroidThermoPhysicalModelResult
    times  ::Vector{Float64}
    E_in   ::Vector{Float64}
    E_out  ::Vector{Float64}
    force  ::Vector{SVector{3, Float64}}
    torque ::Vector{SVector{3, Float64}}

    times_to_save          ::Vector{Float64}
    depth_nodes            ::Vector{Float64}
    surface_temperature    ::Matrix{Float64}
    subsurface_temperature ::Dict{Int, Matrix{Float64}}
    face_forces            ::Matrix{SVector{3, Float64}}
end


"""
Outer constructor of `SingleAsteroidThermoPhysicalModelResult`

# Arguments
- `stpm`          : Thermophysical model for a single asteroid
- `ephem`         : Ephemerides
- `times_to_save` : Timesteps to save temperature
- `face_ID`       : Face indices to save subsurface temperature
"""
function SingleAsteroidThermoPhysicalModelResult(stpm::SingleAsteroidThermoPhysicalModel, ephem, times_to_save::Vector{Float64}, face_ID::Vector{Int})
    n_step = length(ephem.time)             # Number of time steps
    n_step_to_save = length(times_to_save)  # Number of time steps to save temperature
    n_face = length(stpm.shape.faces)       # Number of faces of the shape model

    E_in   = zeros(n_step)
    E_out  = zeros(n_step)
    force  = zeros(SVector{3, Float64}, n_step)
    torque = zeros(SVector{3, Float64}, n_step)

    depth_nodes = stpm.thermo_params.Δz * (0:stpm.thermo_params.n_depth-1)
    surface_temperature = zeros(n_face, n_step_to_save)
    subsurface_temperature = Dict{Int,Matrix{Float64}}(
        i => zeros(stpm.thermo_params.n_depth, n_step_to_save) for i in face_ID
    )
    face_forces = zeros(SVector{3, Float64}, n_face, n_step_to_save)

    return SingleAsteroidThermoPhysicalModelResult(
        ephem.time,
        E_in,
        E_out,
        force,
        torque,
        times_to_save,
        depth_nodes,
        surface_temperature,
        subsurface_temperature,
        face_forces,
    )
end


"""
    struct BinaryAsteroidThermoPhysicalModelResult

Output data format for `BinaryAsteroidThermoPhysicalModel`

# Fields
- `pri` : TPM result for the primary
- `sec` : TPM result for the secondary
"""
struct BinaryAsteroidThermoPhysicalModelResult
    pri::SingleAsteroidThermoPhysicalModelResult
    sec::SingleAsteroidThermoPhysicalModelResult
end


"""
Outer constructor of `BinaryAsteroidThermoPhysicalModelResult`

# Arguments
- `btpm`          : Thermophysical model for a binary asteroid
- `ephem`         : Ephemerides
- `times_to_save` : Timesteps to save temperature (Common to both the primary and the secondary)
- `face_ID_pri`   : Face indices to save subsurface temperature of the primary
- `face_ID_sec`   : Face indices to save subsurface temperature of the secondary
"""
function BinaryAsteroidThermoPhysicalModelResult(btpm::BinaryAsteroidThermoPhysicalModel, ephem, times_to_save::Vector{Float64}, face_ID_pri::Vector{Int}, face_ID_sec::Vector{Int})
    result_pri = SingleAsteroidThermoPhysicalModelResult(btpm.pri, ephem, times_to_save, face_ID_pri)
    result_sec = SingleAsteroidThermoPhysicalModelResult(btpm.sec, ephem, times_to_save, face_ID_sec)

    return BinaryAsteroidThermoPhysicalModelResult(result_pri, result_sec)
end


"""
    update_TPM_result!(result::SingleAsteroidThermoPhysicalModelResult, stpm::SingleAsteroidThermoPhysicalModel, i_time::Integer)

Save the results of TPM at the time step `i_time` to `result`.

# Arguments
- `result` : Output data format for `SingleAsteroidThermoPhysicalModel`
- `stpm`   : Thermophysical model for a single asteroid
- `i_time` : Time step to save data
"""
function update_TPM_result!(result::SingleAsteroidThermoPhysicalModelResult, stpm::SingleAsteroidThermoPhysicalModel, i_time::Integer)
    result.E_in[i_time]   = energy_in(stpm)
    result.E_out[i_time]  = energy_out(stpm)
    result.force[i_time]  = stpm.force
    result.torque[i_time] = stpm.torque

    t  = result.times[i_time]  # Current time

    if t in result.times_to_save  # In the step of saving temperature
        i_time_save = findfirst(isequal(t), result.times_to_save)

        result.surface_temperature[:, i_time_save] .= surface_temperature(stpm)

        for (i, temperature) in result.subsurface_temperature
            temperature[:, i_time_save] .= stpm.temperature[:, i]
        end

        result.face_forces[:, i_time_save] .= stpm.face_forces
    end
end


"""
    update_TPM_result!(result::BinaryAsteroidThermoPhysicalModelResult, btpm::BinaryAsteroidThermoPhysicalModel, ephem, i_time::Integer)

Save the results of TPM at the time step `i_time` to `result`.

# Arguments
- `result` : Output data format for `BinaryAsteroidThermoPhysicalModel`
- `btpm`   : Thermophysical model for a binary asteroid
- `ephem`  : Ephemerides
- `i_time`     : Time step
"""
function update_TPM_result!(result::BinaryAsteroidThermoPhysicalModelResult, btpm::BinaryAsteroidThermoPhysicalModel, i_time::Integer)
    update_TPM_result!(result.pri, btpm.pri, i_time)
    update_TPM_result!(result.sec, btpm.sec, i_time)
end


"""
    export_TPM_results(dirpath, result::SingleAsteroidThermoPhysicalModelResult)

Export the result of `SingleAsteroidThermoPhysicalModel` to CSV files. 
The output files are saved in the following directory structure:

    dirpath
    ├── physical_quantities.csv
    ├── subsurface_temperature.csv
    ├── surface_temperature.csv
    └── thermal_force.csv

# Arguments
- `dirpath` : Path to the directory to save CSV files.
- `result`  : Output data format for `SingleAsteroidThermoPhysicalModel`
"""
function export_TPM_results(dirpath, result::SingleAsteroidThermoPhysicalModelResult)
    
    df = DataFrame()
    df.time     = result.times
    df.E_in     = result.E_in
    df.E_out    = result.E_out
    df.force_x  = [f[1] for f in result.force]
    df.force_y  = [f[2] for f in result.force]
    df.force_z  = [f[3] for f in result.force]
    df.torque_x = [τ[1] for τ in result.torque]
    df.torque_y = [τ[2] for τ in result.torque]
    df.torque_z = [τ[3] for τ in result.torque]
    
    CSV.write(joinpath(dirpath, "physical_quantities.csv"), df)

    ##= Surface temperature =##
    filepath = joinpath(dirpath, "surface_temperature.csv")
    df = hcat(
        DataFrame(time = result.times_to_save),
        DataFrame(
            result.surface_temperature',
            ["face_$(i)" for i in 1:size(result.surface_temperature, 1)]
        ),
    )

    CSV.write(filepath, df)

    ##= Subsurface temperature =##
    filepath = joinpath(dirpath, "subsurface_temperature.csv")
    
    nrows = length(result.depth_nodes) * length(result.times_to_save)
    df = DataFrame(
        time  = reshape([t for _ in result.depth_nodes, t in result.times_to_save], nrows),
        depth = reshape([d for d in result.depth_nodes, _ in result.times_to_save], nrows),
    )

    # Add a column for each face
    for (i, subsurface_temperature) in collect(result.subsurface_temperature)
        df[:, "face_$(i)"] =
            reshape(subsurface_temperature, length(subsurface_temperature))
    end

    # Sort the columns by the face ID
    keys_sorted = sort(names(df[:, 3:end]), by=x->parse(Int, replace(x, r"[^0-9]" => "")))
    df = df[:, ["time", "depth", keys_sorted...]]

    CSV.write(filepath, df)

    ##= Thermal force on every face of the shape model =##
    filepath = joinpath(dirpath, "thermal_force.csv")

    n_face = size(result.face_forces, 1)  # Number of faces of the shape model
    n_step = size(result.face_forces, 2)  # Number of time steps to save temperature
    nrows = n_face * n_step

    df = DataFrame(
        time = reshape([t for _ in 1:n_face, t in result.times_to_save], nrows),
        face = reshape([i for i in 1:n_face, _ in result.times_to_save], nrows),
    )
    df.x = reshape([f[1] for f in result.face_forces], nrows)  # x-component of the thermal force
    df.y = reshape([f[2] for f in result.face_forces], nrows)  # y-component of the thermal force
    df.z = reshape([f[3] for f in result.face_forces], nrows)  # z-component of the thermal force

    CSV.write(filepath, df)
end


"""
    export_TPM_results(dirpath, result::BinaryAsteroidThermoPhysicalModelResult)

Export the result of `BinaryAsteroidThermoPhysicalModel` to CSV files. 
The output files are saved in the following directory structure:

    dirpath
    ├── pri
    │   ├── physical_quantities.csv
    │   ├── subsurface_temperature.csv
    │   ├── surface_temperature.csv
    │   └── thermal_force.csv
    └── sec
        ├── physical_quantities.csv
        ├── subsurface_temperature.csv
        ├── surface_temperature.csv
        └── thermal_force.csv

# Arguments
- `dirpath` : Path to the directory to save CSV files.
- `result`  : Output data format for `BinaryAsteroidThermoPhysicalModel`
"""
function export_TPM_results(dirpath, result::BinaryAsteroidThermoPhysicalModelResult)
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
    init_temperature!(btpm::BinaryTBinaryAsteroidThermoPhysicalModelPM, T₀::Real)

Initialize all temperature cells at the given temperature `T₀`

# Arguments
- `btpm` : Thermophysical model for a binary asteroid
- `T₀`   : Initial temperature of all cells [K]
"""
function init_temperature!(btpm::BinaryAsteroidThermoPhysicalModel, T₀::Real)
    init_temperature!(btpm.pri, T₀)
    init_temperature!(btpm.sec, T₀)
end


# ****************************************************************
#                    Thermophysical modeling
# ****************************************************************


"""
    run_TPM!(stpm::SingleAsteroidThermoPhysicalModel, ephem, savepath)

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
function run_TPM!(stpm::SingleAsteroidThermoPhysicalModel, ephem, times_to_save::Vector{Float64}, face_ID::Vector{Int}; show_progress=true)

    result = SingleAsteroidThermoPhysicalModelResult(stpm, ephem, times_to_save, face_ID)

    ## ProgressMeter setting
    if show_progress
        p = Progress(length(ephem.time); dt=1, desc="Running TPM...", showspeed=true)
        ProgressMeter.ijulia_behavior(:clear)
    end
    
    for i_time in eachindex(ephem.time)
        r☉ = ephem.sun[i_time]

        update_flux_sun!(stpm, r☉)
        update_flux_scat_single!(stpm)
        update_flux_rad_single!(stpm)
        
        update_thermal_force!(stpm)

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
    run_TPM!(btpm::BinaryAsteroidThermoPhysicalModel, ephem, savepath)

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
function run_TPM!(btpm::BinaryAsteroidThermoPhysicalModel, ephem, times_to_save::Vector{Float64}, face_ID_pri::Vector{Int}, face_ID_sec::Vector{Int}; show_progress=true)

    result = BinaryAsteroidThermoPhysicalModelResult(btpm, ephem, times_to_save, face_ID_pri, face_ID_sec)

    ## ProgressMeter setting
    if show_progress
        p = Progress(length(ephem.time); dt=1, desc="Running TPM...", showspeed=true)
        ProgressMeter.ijulia_behavior(:clear)
    end
    
    for i_time in eachindex(ephem.time)
        r☉₁ = ephem.sun1[i_time]  # Sun's position in the primary's frame
        r☉₂ = ephem.sun2[i_time]  # Sun's position in the secondary's frame
        rₛ  = ephem.sec[i_time]   # Secondary's position in the primary's frame
        R₂₁ = ephem.S2P[i_time]   # Rotation matrix from secondary to primary frames

        ## Update enegey flux
        update_flux_sun!(btpm, r☉₁, r☉₂)
        mutual_shadowing!(btpm, r☉₁, rₛ, R₂₁)  # Mutual-shadowing (eclipse)
        update_flux_scat_single!(btpm)
        update_flux_rad_single!(btpm)
        mutual_heating!(btpm, rₛ, R₂₁)         # Mutual-heating

        update_thermal_force!(btpm)

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
