#=
tpm_state.jl

Type definitions for single and binary asteroid thermophysical simulation state.
=#


"""
Abstract type for asteroid thermophysical simulation state.
"""
abstract type AbstractAsteroidThermoPhysicalState end


"""
    struct SingleAsteroidThermoPhysicalState <: AbstractAsteroidThermoPhysicalState

Internal simulation state for a single-asteroid thermophysical model.
Holds the mutable arrays that evolve during a `solve` call.
The problem definition (shape, parameters, flags, boundary conditions) is
accessed via the `problem` field to avoid duplication.

# Fields
- `problem`           : Problem definition (shape, thermo_params, flags, BCs)
- `solver_cache`      : Pre-allocated cache for the heat-conduction solver
- `illuminated_faces` : Illumination flag for each face
- `flux_sun`          : Direct solar flux on each face [W/m²]
- `flux_scat`         : Scattered-light flux on each face [W/m²]
- `flux_rad`          : Thermal-emission flux from surrounding faces [W/m²]
- `temperature`       : Temperature matrix `(n_depth, n_face)` [K]
- `face_forces`       : Thermal recoil force on each face [N]
- `force`             : Net thermal recoil force in body-fixed frame [N]
- `torque`            : Net thermal recoil torque in body-fixed frame [N⋅m]
"""
struct SingleAsteroidThermoPhysicalState{
    Pr  <: SingleAsteroidThermoPhysicalProblem,
    HCC <: HeatConductionCache,
} <: AbstractAsteroidThermoPhysicalState
    problem           ::Pr
    solver_cache      ::HCC
    illuminated_faces ::Vector{Bool}
    flux_sun          ::Vector{Float64}
    flux_scat         ::Vector{Float64}
    flux_rad          ::Vector{Float64}
    temperature       ::Matrix{Float64}  # (n_depth, n_face)
    face_forces       ::Vector{SVector{3, Float64}}
    force             ::MVector{3, Float64}
    torque            ::MVector{3, Float64}
end


"""
    struct BinaryAsteroidThermoPhysicalState <: AbstractAsteroidThermoPhysicalState

Internal simulation state for a binary-asteroid thermophysical model.

# Fields
- `problem`   : Binary problem definition (mutual shadowing/heating flags)
- `primary`   : Simulation state for the primary body
- `secondary` : Simulation state for the secondary body

# Invariant
The inner constructor enforces:
```
state.primary.problem   === state.problem.primary
state.secondary.problem === state.problem.secondary
```
Use `_build_binary_state` rather than constructing directly to guarantee consistency.
"""
struct BinaryAsteroidThermoPhysicalState{
    St1 <: SingleAsteroidThermoPhysicalState,
    St2 <: SingleAsteroidThermoPhysicalState,
} <: AbstractAsteroidThermoPhysicalState
    problem   ::BinaryAsteroidThermoPhysicalProblem
    primary   ::St1
    secondary ::St2

    function BinaryAsteroidThermoPhysicalState(
        problem   ::BinaryAsteroidThermoPhysicalProblem,
        primary   ::St1,
        secondary ::St2,
    ) where {St1 <: SingleAsteroidThermoPhysicalState, St2 <: SingleAsteroidThermoPhysicalState}
        primary.problem   === problem.primary   || error("primary.problem ≢ problem.primary: use _build_binary_state to construct")
        secondary.problem === problem.secondary || error("secondary.problem ≢ problem.secondary: use _build_binary_state to construct")
        new{St1, St2}(problem, primary, secondary)
    end
end


"""
    struct HierarchicalSingleAsteroidThermoPhysicalState <: AbstractAsteroidThermoPhysicalState

Internal simulation state for a single-asteroid thermophysical model with surface roughness.
Mirrors the two-level structure of `HierarchicalShapeModel`: a global-level state identical in
layout to `SingleAsteroidThermoPhysicalState`, plus per-face sub-states for the roughness models.

# Fields
## Global level
- `problem`           : Problem definition (shape must be `HierarchicalShapeModel`)
- `solver_cache`      : Pre-allocated cache for the global heat-conduction solver
- `illuminated_faces` : Illumination flag for each global face
- `flux_sun`          : Direct solar flux on each global face [W/m²]
- `flux_scat`         : Scattered-light flux on each global face [W/m²]
- `flux_rad`          : Thermal-emission flux from surrounding global faces [W/m²]
- `temperature`       : Temperature matrix `(n_depth, n_global_faces)` [K]
- `face_forces`       : Thermal recoil force on each global face [N]
- `force`             : Net thermal recoil force in body-fixed frame [N]
- `torque`            : Net thermal recoil torque in body-fixed frame [N⋅m]
## Sub-face level (symmetric with `HierarchicalShapeModel`)
- `face_roughness_indices` : Maps global face index → `roughness_states` index (0 = no roughness);
                             length = `n_global_faces`; same direction as `HierarchicalShapeModel.face_roughness_indices`
- `roughness_states`       : Independent `SingleAsteroidThermoPhysicalState` per roughness-carrying face;
                             length = number of global faces that have a roughness model
"""
struct HierarchicalSingleAsteroidThermoPhysicalState{
    Pr  <: SingleAsteroidThermoPhysicalProblem{<:HierarchicalShapeModel},
    HCC <: HeatConductionCache,
    St  <: SingleAsteroidThermoPhysicalState,
} <: AbstractAsteroidThermoPhysicalState
    problem           ::Pr
    solver_cache      ::HCC
    illuminated_faces ::Vector{Bool}
    flux_sun          ::Vector{Float64}
    flux_scat         ::Vector{Float64}
    flux_rad          ::Vector{Float64}
    temperature       ::Matrix{Float64}  # (n_depth, n_global_faces)
    face_forces       ::Vector{SVector{3, Float64}}
    force             ::MVector{3, Float64}
    torque            ::MVector{3, Float64}
    face_roughness_indices ::Vector{Int}
    roughness_states       ::Vector{St}
end


"""
    surface_temperature(state::SingleAsteroidThermoPhysicalState) -> T_surface

Extract the surface temperature (uppermost layer) for all faces.

# Returns
- `T_surface::Vector{Float64}` : Surface temperature for each face [K]
"""
surface_temperature(state::SingleAsteroidThermoPhysicalState) = state.temperature[begin, :]


"""
    surface_temperature(state::HierarchicalSingleAsteroidThermoPhysicalState) -> T_surface

Extract the global-level surface temperature (uppermost layer) for all global faces.

# Returns
- `T_surface::Vector{Float64}` : Surface temperature for each global face [K]
"""
surface_temperature(state::HierarchicalSingleAsteroidThermoPhysicalState) = state.temperature[begin, :]
