#=
TPM.jl

Core types and functions for thermophysical modeling of asteroids.
This file defines the main data structures and simulation functions including:
- Heat conduction solver types (Explicit/Implicit Euler, Crank-Nicolson)
- Boundary condition types (Radiation, Insulation, Isothermal)
- Single and binary asteroid thermophysical model structures
- Temperature initialization and evolution functions
- Result storage and export functionality
=#

# ╔═══════════════════════════════════════════════════════════════════╗
# ║          Types of solvers for heat conduction equations           ║
# ╚═══════════════════════════════════════════════════════════════════╝

"""
Abstract type of a solver for a heat conduction equation 
"""
abstract type HeatConductionSolver end


"""
Type of the explicit (forward) Euler method:
- Explicit in time
- First order in time

The `ExplicitEulerSolver` type includes a vector for the temperature at the next time step.
"""
struct ExplicitEulerSolver <: HeatConductionSolver
    x::Vector{Float64}  # Temperature vector for the next time step
end

ExplicitEulerSolver(thermo_params::AbstractThermoParams) = ExplicitEulerSolver(thermo_params.n_depth)
ExplicitEulerSolver(N::Integer) = ExplicitEulerSolver(zeros(N))


"""
Type of the implicit (backward) Euler method:
- Implicit in time (Unconditionally stable in the heat conduction equation)
- First order in time

The `ImplicitEulerSolver` type has vectors for the tridiagonal matrix algorithm.
"""
struct ImplicitEulerSolver <: HeatConductionSolver
    a::Vector{Float64}
    b::Vector{Float64}
    c::Vector{Float64}
    d::Vector{Float64}
    x::Vector{Float64}
end

ImplicitEulerSolver(thermo_params::AbstractThermoParams) = ImplicitEulerSolver(thermo_params.n_depth)
ImplicitEulerSolver(N::Integer) = ImplicitEulerSolver(zeros(N), zeros(N), zeros(N), zeros(N), zeros(N))


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


# ╔═══════════════════════════════════════════════════════════════════╗
# ║    Types of boundary conditions for heat conduction equations    ║
# ╚═══════════════════════════════════════════════════════════════════╝

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


# ╔═══════════════════════════════════════════════════════════════════╗
# ║                Types of thermophsycal models                      ║
# ╚═══════════════════════════════════════════════════════════════════╝


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
- `shape`         : Shape model
- `thermo_params` : Thermophysical parameters

- `illuminated_faces` : Boolean array indicating illumination state for each face (used for batch processing with AsteroidShapeModels.jl v0.4.0+)
- `flux_sun`          : Flux of direct sunlight on each face [W/m²]
- `flux_scat`         : Flux of scattered light on each face [W/m²]
- `flux_rad`          : Flux of thermal emission from surrounding surface on each face [W/m²]
- `temperature`       : Temperature matrix `(n_depth, n_face)` according to the number of depth cells `n_depth` and the number of faces `n_face`.

- `face_forces` : Thermal force on each face
- `force`       : Thermal recoil force at body-fixed frame (Yarkovsky effect)
- `torque`      : Thermal recoil torque at body-fixed frame (YORP effect)

- `SELF_SHADOWING` : Flag to consider self-shadowing
- `SELF_HEATING`   : Flag to consider self-heating
- `SOLVER`         : Solver of heat conduction equation
- `BC_UPPER`       : Boundary condition at the upper boundary
- `BC_LOWER`       : Boundary condition at the lower boundary

# TODO:
- roughness_maps   ::ShapeModel[]
"""
struct SingleAsteroidThermoPhysicalModel{P<:AbstractThermoParams, S<:HeatConductionSolver, BU<:BoundaryCondition, BL<:BoundaryCondition} <: AbstractAsteroidThermoPhysicalModel
    shape         ::ShapeModel
    thermo_params ::P

    illuminated_faces ::Vector{Bool}     # Illumination state for each face (for batch processing)
    flux_sun          ::Vector{Float64}  # Flux of direct sunlight
    flux_scat         ::Vector{Float64}  # Flux of scattered light
    flux_rad          ::Vector{Float64}  # Flux of thermal emission from surrounding surface
    temperature       ::Matrix{Float64}  # (n_depth, n_face)

    face_forces ::Vector{SVector{3, Float64}}
    force       ::MVector{3, Float64}
    torque      ::MVector{3, Float64}

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

# Notes
- If `SELF_SHADOWING` is true and face_visibility_graph is not built, it will be automatically built
- This may take some time for large shape models
- To avoid automatic building, pre-build with `build_face_visibility_graph!` or load shape with `with_face_visibility=true`
"""
function SingleAsteroidThermoPhysicalModel(shape, thermo_params; SELF_SHADOWING, SELF_HEATING, SOLVER, BC_UPPER, BC_LOWER)
    # Automatically build face_visibility_graph if needed for self-shadowing
    if SELF_SHADOWING
        if isnothing(shape.face_visibility_graph)
            @info "Building face_visibility_graph for self-shadowing..."
            build_face_visibility_graph!(shape)
        end
    end

    broadcast_thermo_params!(thermo_params, shape)

    n_depth = thermo_params.n_depth
    n_face = length(shape.faces)

    illuminated_faces = zeros(Bool, n_face)
    flux_sun = zeros(n_face)
    flux_scat = zeros(n_face)
    flux_rad = zeros(n_face)
    temperature = zeros(n_depth, n_face)

    face_forces = zeros(SVector{3, Float64}, n_face)
    force  = zero(MVector{3, Float64})
    torque = zero(MVector{3, Float64})

    SingleAsteroidThermoPhysicalModel(shape, thermo_params, illuminated_faces, flux_sun, flux_scat, flux_rad, temperature, face_forces, force, torque, SELF_SHADOWING, SELF_HEATING, SOLVER, BC_UPPER, BC_LOWER)
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
    BinaryAsteroidThermoPhysicalModel(pri, sec; MUTUAL_SHADOWING, MUTUAL_HEATING) -> btpm

Construct a thermophysical model for a binary asteroid (`BinaryAsteroidThermoPhysicalModel`).

# Arguments
- `pri` : TPM for the primary asteroid
- `sec` : TPM for the secondary asteroid

# Keyword Arguments
- `MUTUAL_SHADOWING` : Flag to consider mutual shadowing (required)
- `MUTUAL_HEATING`   : Flag to consider mutual heating (required)

# Notes
- Both `MUTUAL_SHADOWING` and `MUTUAL_HEATING` must be explicitly specified
- If `MUTUAL_SHADOWING` is true and BVH is not built, it will be automatically built
- This may take some time for large shape models
- To avoid automatic BVH building, pre-build with `build_bvh!` or load shapes with `with_bvh=true`

# Example
```julia
# Explicit specification required
btpm = BinaryAsteroidThermoPhysicalModel(pri, sec; 
    MUTUAL_SHADOWING = true,
    MUTUAL_HEATING   = true,
)

# Or disable mutual effects
btpm = BinaryAsteroidThermoPhysicalModel(pri, sec; 
    MUTUAL_SHADOWING = false,
    MUTUAL_HEATING   = false,
)
```
"""
function BinaryAsteroidThermoPhysicalModel(pri, sec; MUTUAL_SHADOWING, MUTUAL_HEATING)
    # Automatically build BVH if needed for mutual shadowing
    if MUTUAL_SHADOWING
        if isnothing(pri.shape.bvh)
            @info "Building BVH for primary shape..."
            build_bvh!(pri.shape)
        end
        if isnothing(sec.shape.bvh)
            @info "Building BVH for secondary shape..."
            build_bvh!(sec.shape)
        end
    end
    
    BinaryAsteroidThermoPhysicalModel(pri, sec, MUTUAL_SHADOWING, MUTUAL_HEATING)
end


"""
    surface_temperature(stpm::SingleAsteroidThermoPhysicalModel) -> T_surface

Extract the surface temperature values for all faces of the asteroid.

# Arguments
- `stpm::SingleAsteroidThermoPhysicalModel` : Thermophysical model for a single asteroid

# Returns
- `T_surface::Vector{Float64}` : Surface temperature for each face [K]

# Notes
- Returns temperatures from the uppermost layer (index 1) of the temperature matrix
- The length of the returned vector equals the number of faces in the shape model
- Useful for visualization, thermal emission calculations, and analysis

# Example
```julia
T_surf = surface_temperature(stpm)
T_mean = mean(T_surf)  # Average surface temperature
T_max = maximum(T_surf)  # Hottest point
T_min = minimum(T_surf)  # Coldest point
```

# See Also
- `stpm.temperature` for the full temperature matrix including subsurface
"""
surface_temperature(stpm::SingleAsteroidThermoPhysicalModel) = stpm.temperature[begin, :]


# ╔═══════════════════════════════════════════════════════════════════╗
# ║                     Output data format                            ║
# ╚═══════════════════════════════════════════════════════════════════╝

"""
    struct SingleAsteroidThermoPhysicalModelResult

Output data format for `SingleAsteroidThermoPhysicalModel` 

# Fields
## Saved at all time steps
- `times`  : Timesteps, given the same vector as `ephem.time` [s]
- `E_in`   : Input energy per second on the whole surface [W]
- `E_out`  : Output energy per second from the whole surface [W]
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
