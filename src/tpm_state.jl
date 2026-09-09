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
## Surface roughness (empty for a smooth surface)
- `face_roughness_indices` : Maps face index → `roughness_states` index (0 = no roughness);
                             length = `n_face` when the shape has roughness, empty otherwise.
                             Mirrors `shape.roughness.face_roughness_indices`, but maps to
                             *independent* per-face states rather than shared models.
- `roughness_states`       : Independent sub-face state per roughness-carrying face;
                             empty when the shape has no roughness. Each sub-state is itself a
                             `SingleAsteroidThermoPhysicalState` with empty roughness (the
                             roughness models are smooth `ShapeModel`s)

# Notes
When the shape carries surface roughness (`has_roughness(problem.shape)`), the fields above
describe the global faces, and each global face with a roughness model additionally has its
own sub-face state in `roughness_states`. Code that iterates the roughness states can simply
loop over them: for a smooth surface the vectors are empty and the loop does nothing.
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
    face_roughness_indices ::Vector{Int}
    roughness_states       ::Vector{SingleAsteroidThermoPhysicalState{Pr, HCC}}
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
    surface_temperature(state::SingleAsteroidThermoPhysicalState) -> T_surface

Extract the surface temperature (uppermost layer) for all faces.
For a shape with surface roughness, these are the global faces; the sub-face
temperatures live in `state.roughness_states`.

# Returns
- `T_surface::Vector{Float64}` : Surface temperature for each face [K]
"""
surface_temperature(state::SingleAsteroidThermoPhysicalState) = state.temperature[begin, :]
