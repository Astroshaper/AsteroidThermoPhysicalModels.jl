#=
tpm_init.jl

Temperature initialization for thermophysical models.
=#


"""
    init_temperature!(stpm::SingleAsteroidThermoPhysicalModel, T₀::Real)

Initialize all temperature cells at the uniform temperature `T₀`.

# Arguments
- `stpm` : Internal thermophysical model
- `T₀`   : Initial temperature [K]
"""
function init_temperature!(stpm::SingleAsteroidThermoPhysicalModel, T₀::Real)
    stpm.temperature .= T₀
end


"""
    init_temperature!(stpm::SingleAsteroidThermoPhysicalModel, T₀::AbstractMatrix)

Initialize temperatures from a full depth–face temperature matrix.
The matrix must have size `(n_depth, n_face)`, matching `stpm.temperature`.

# Arguments
- `stpm` : Internal thermophysical model
- `T₀`   : Temperature matrix of size `(n_depth, n_face)` [K]
"""
function init_temperature!(stpm::SingleAsteroidThermoPhysicalModel, T₀::AbstractMatrix)
    stpm.temperature .= T₀
end


"""
    init_temperature!(btpm::BinaryAsteroidThermoPhysicalModel, T₀::Real)

Initialize all temperature cells in both bodies at the uniform temperature `T₀`.

# Arguments
- `btpm` : Thermophysical model for a binary asteroid
- `T₀`   : Initial temperature of all cells [K]
"""
function init_temperature!(btpm::BinaryAsteroidThermoPhysicalModel, T₀::Real)
    init_temperature!(btpm.pri, T₀)
    init_temperature!(btpm.sec, T₀)
end


"""
    init_temperature!(btpm::BinaryAsteroidThermoPhysicalModel, T₀_primary, T₀_secondary)

Initialize temperatures with separate values for the primary and secondary bodies.

Each argument can be a `Real` (uniform) or an `AbstractMatrix` of size `(n_depth, n_face)`.

# Arguments
- `btpm`         : Thermophysical model for a binary asteroid
- `T₀_primary`   : Initial temperature for the primary body [K]
- `T₀_secondary` : Initial temperature for the secondary body [K]
"""
function init_temperature!(btpm::BinaryAsteroidThermoPhysicalModel, T₀_primary, T₀_secondary)
    init_temperature!(btpm.pri, T₀_primary)
    init_temperature!(btpm.sec, T₀_secondary)
end
