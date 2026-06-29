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
    Pr <: SingleAsteroidThermoPhysicalProblem,
    S  <: HeatConductionCache,
} <: AbstractAsteroidThermoPhysicalState
    problem           ::Pr
    solver_cache      ::S
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
    S1 <: SingleAsteroidThermoPhysicalState,
    S2 <: SingleAsteroidThermoPhysicalState,
} <: AbstractAsteroidThermoPhysicalState
    problem   ::BinaryAsteroidThermoPhysicalProblem
    primary   ::S1
    secondary ::S2

    function BinaryAsteroidThermoPhysicalState(
        problem   ::BinaryAsteroidThermoPhysicalProblem,
        primary   ::S1,
        secondary ::S2,
    ) where {S1 <: SingleAsteroidThermoPhysicalState, S2 <: SingleAsteroidThermoPhysicalState}
        primary.problem   === problem.primary   || error("primary.problem ≢ problem.primary: use _build_binary_state to construct")
        secondary.problem === problem.secondary || error("secondary.problem ≢ problem.secondary: use _build_binary_state to construct")
        new{S1, S2}(problem, primary, secondary)
    end
end


"""
    surface_temperature(state::SingleAsteroidThermoPhysicalState) -> T_surface

Extract the surface temperature (uppermost layer) for all faces.

# Returns
- `T_surface::Vector{Float64}` : Surface temperature for each face [K]
"""
surface_temperature(state::SingleAsteroidThermoPhysicalState) = state.temperature[begin, :]
