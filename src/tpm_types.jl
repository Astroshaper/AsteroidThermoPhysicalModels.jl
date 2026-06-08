#=
tpm_types.jl

Type definitions for single and binary asteroid thermophysical models.
=#

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
