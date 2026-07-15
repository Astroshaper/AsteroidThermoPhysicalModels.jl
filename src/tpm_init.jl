#=
tpm_init.jl

Temperature initialization for thermophysical simulation states.
=#


"""
    init_temperature!(state::SingleAsteroidThermoPhysicalState, T₀::Real)

Initialize all temperature cells at the uniform temperature `T₀`.

# Arguments
- `state` : Thermophysical simulation state for a single asteroid
- `T₀`   : Initial temperature [K]
"""
function init_temperature!(state::SingleAsteroidThermoPhysicalState, T₀::Real)
    state.temperature .= T₀
end


"""
    init_temperature!(state::SingleAsteroidThermoPhysicalState, T₀::AbstractMatrix)

Initialize temperatures from a full depth–face temperature matrix.
The matrix must have size `(n_depth, n_face)`, matching `state.temperature`.

# Arguments
- `state` : Thermophysical simulation state for a single asteroid
- `T₀`   : Temperature matrix of size `(n_depth, n_face)` [K]
"""
function init_temperature!(state::SingleAsteroidThermoPhysicalState, T₀::AbstractMatrix)
    state.temperature .= T₀
end


"""
    init_temperature!(state::HierarchicalSingleAsteroidThermoPhysicalState, T₀::Real)

Initialize all temperature cells at the uniform temperature `T₀`, including all sub-face states.

# Arguments
- `state` : Thermophysical simulation state for a single asteroid with surface roughness
- `T₀`   : Initial temperature [K]
"""
function init_temperature!(state::HierarchicalSingleAsteroidThermoPhysicalState, T₀::Real)
    state.temperature .= T₀
    for rs in state.roughness_states
        rs.temperature .= T₀
    end
end


"""
    init_temperature!(state::HierarchicalSingleAsteroidThermoPhysicalState, T₀::AbstractMatrix)

Initialize global temperatures from a depth–face matrix.
Each sub-face state is initialized to the surface temperature of its parent global face.

# Arguments
- `state` : Thermophysical simulation state for a single asteroid with surface roughness
- `T₀`   : Temperature matrix of size `(n_depth, n_global_faces)` [K]
"""
function init_temperature!(state::HierarchicalSingleAsteroidThermoPhysicalState, T₀::AbstractMatrix)
    state.temperature .= T₀
    for (i, k) in enumerate(state.face_roughness_indices)
        k == 0 && continue
        state.roughness_states[k].temperature .= T₀[begin, i]
    end
end


"""
    init_temperature!(state::BinaryAsteroidThermoPhysicalState, T₀::Real)

Initialize all temperature cells in both bodies at the uniform temperature `T₀`.

# Arguments
- `state` : Thermophysical simulation state for a binary asteroid
- `T₀`   : Initial temperature of all cells [K]
"""
function init_temperature!(state::BinaryAsteroidThermoPhysicalState, T₀::Real)
    init_temperature!(state.primary, T₀)
    init_temperature!(state.secondary, T₀)
end


"""
    init_temperature!(state::BinaryAsteroidThermoPhysicalState, T₀_primary, T₀_secondary)

Initialize temperatures with separate values for the primary and secondary bodies.

Each argument can be a `Real` (uniform) or an `AbstractMatrix` of size `(n_depth, n_face)`.

# Arguments
- `state`         : Thermophysical simulation state for a binary asteroid
- `T₀_primary`   : Initial temperature for the primary body [K]
- `T₀_secondary` : Initial temperature for the secondary body [K]
"""
function init_temperature!(state::BinaryAsteroidThermoPhysicalState, T₀_primary, T₀_secondary)
    init_temperature!(state.primary, T₀_primary)
    init_temperature!(state.secondary, T₀_secondary)
end
