#=
ephem_types.jl

Formal ephemerides types for thermophysical simulations.
These replace the informal NamedTuple previously used as `ephem`.

Type hierarchy:
    AbstractAsteroidEphemerides
    ├── AbstractSingleAsteroidEphemerides
    │   └── SingleAsteroidEphemerides{R}
    └── AbstractBinaryAsteroidEphemerides
        └── BinaryAsteroidEphemerides{R}

Type parameter R:
    R = Nothing                          : temperature only (no force/torque)
    R = Vector{SMatrix{3,3,Float64,9}}   : force/torque output in the inertial frame
=#

# ╔═══════════════════════════════════════════════════════════════════╗
# ║                       Abstract types                              ║
# ╚═══════════════════════════════════════════════════════════════════╝

"""
    abstract type AbstractAsteroidEphemerides

Supertype for all ephemerides types used in thermophysical simulations.
"""
abstract type AbstractAsteroidEphemerides end

"""
    abstract type AbstractSingleAsteroidEphemerides <: AbstractAsteroidEphemerides

Supertype for single-asteroid ephemerides.
"""
abstract type AbstractSingleAsteroidEphemerides <: AbstractAsteroidEphemerides end

"""
    abstract type AbstractBinaryAsteroidEphemerides <: AbstractAsteroidEphemerides

Supertype for binary-asteroid ephemerides.
"""
abstract type AbstractBinaryAsteroidEphemerides <: AbstractAsteroidEphemerides end

"""
    Base.length(ephem::AbstractAsteroidEphemerides) -> Int

Return the number of timesteps in the ephemerides.
"""
Base.length(ephem::AbstractAsteroidEphemerides) = length(ephem.times)

# ╔═══════════════════════════════════════════════════════════════════╗
# ║                     Single asteroid                               ║
# ╚═══════════════════════════════════════════════════════════════════╝

"""
    struct SingleAsteroidEphemerides{R} <: AbstractSingleAsteroidEphemerides

Ephemerides for a single-asteroid thermophysical simulation.

The type parameter `R` controls whether force and torque are computed:
- `R = Nothing`                        : temperature only; `R_body_to_inertial = nothing`
- `R = Vector{SMatrix{3,3,Float64,9}}` : force and torque are computed in the inertial frame

# Fields
- `times`              : Simulation timesteps [s]
- `r_sun`              : Sun position vector in the body-fixed frame at each timestep [m]
- `R_body_to_inertial` : Passive rotation matrices from body-fixed to inertial frame,
                         or `nothing` when force/torque are not needed.

# Constructors
    SingleAsteroidEphemerides(times, r_sun)
        -> SingleAsteroidEphemerides{Nothing}
    SingleAsteroidEphemerides(times, r_sun, R_body_to_inertial)
        -> SingleAsteroidEphemerides{typeof(R_body_to_inertial)}
"""
struct SingleAsteroidEphemerides{R} <: AbstractSingleAsteroidEphemerides
    times              ::Vector{Float64}
    r_sun              ::Vector{SVector{3, Float64}}
    R_body_to_inertial ::R

    function SingleAsteroidEphemerides{R}(
        times              ::AbstractVector,
        r_sun              ::AbstractVector,
        R_body_to_inertial ::R,
    ) where {R}
        n = length(times)
        length(r_sun) == n || throw(DimensionMismatch(
            "r_sun ($(length(r_sun))) and times ($n) must have the same length"
        ))
        R_body_to_inertial === nothing || length(R_body_to_inertial) == n || throw(DimensionMismatch(
            "R_body_to_inertial ($(length(R_body_to_inertial))) and times ($n) must have the same length"
        ))
        new{R}(times, r_sun, R_body_to_inertial)
    end
end

# Temperature only: R_body_to_inertial = nothing
SingleAsteroidEphemerides(times, r_sun) =
    SingleAsteroidEphemerides{Nothing}(times, r_sun, nothing)

# With rotation matrices for force/torque computation
SingleAsteroidEphemerides(times, r_sun, R_body_to_inertial) =
    SingleAsteroidEphemerides{typeof(R_body_to_inertial)}(times, r_sun, R_body_to_inertial)


# ╔═══════════════════════════════════════════════════════════════════╗
# ║                     Binary asteroid                               ║
# ╚═══════════════════════════════════════════════════════════════════╝

"""
    struct BinaryAsteroidEphemerides{R} <: AbstractBinaryAsteroidEphemerides

Ephemerides for a binary-asteroid thermophysical simulation.

The type parameter `R` controls whether force and torque are computed:
- `R = Nothing`                        : temperature only; `R_primary_to_inertial = nothing`
- `R = Vector{SMatrix{3,3,Float64,9}}` : force and torque are computed in the inertial frame

The secondary-to-inertial rotation is not stored but can be derived as:

```math
R_{s2i} = R_{p2i} \\cdot R_{p2s}^{\\top}
```

# Fields
- `times`                  : Simulation timesteps [s]
- `r_sun`                  : Sun position vector in the primary body-fixed frame at each timestep [m]
- `r_secondary`            : Secondary position vector in the primary body-fixed frame at each timestep [m]
- `R_primary_to_secondary` : Passive rotation matrices from primary body-fixed to secondary body-fixed frame
- `R_primary_to_inertial`  : Passive rotation matrices from primary body-fixed to inertial frame,
                             or `nothing` when force/torque are not needed.

# Constructors
    BinaryAsteroidEphemerides(times, r_sun, r_secondary, R_primary_to_secondary)
        -> BinaryAsteroidEphemerides{Nothing}
    BinaryAsteroidEphemerides(times, r_sun, r_secondary, R_primary_to_secondary, R_primary_to_inertial)
        -> BinaryAsteroidEphemerides{typeof(R_primary_to_inertial)}
"""
struct BinaryAsteroidEphemerides{R} <: AbstractBinaryAsteroidEphemerides
    times                  ::Vector{Float64}
    r_sun                  ::Vector{SVector{3, Float64}}
    r_secondary            ::Vector{SVector{3, Float64}}
    R_primary_to_secondary ::Vector{SMatrix{3,3,Float64,9}}
    R_primary_to_inertial  ::R

    function BinaryAsteroidEphemerides{R}(
        times                  ::AbstractVector,
        r_sun                  ::AbstractVector,
        r_secondary            ::AbstractVector,
        R_primary_to_secondary ::AbstractVector,
        R_primary_to_inertial  ::R,
    ) where {R}
        n = length(times)
        length(r_sun) == n || throw(DimensionMismatch(
            "r_sun ($(length(r_sun))) and times ($n) must have the same length"
        ))
        length(r_secondary) == n || throw(DimensionMismatch(
            "r_secondary ($(length(r_secondary))) and times ($n) must have the same length"
        ))
        length(R_primary_to_secondary) == n || throw(DimensionMismatch(
            "R_primary_to_secondary ($(length(R_primary_to_secondary))) and times ($n) must have the same length"
        ))
        R_primary_to_inertial === nothing || length(R_primary_to_inertial) == n || throw(DimensionMismatch(
            "R_primary_to_inertial ($(length(R_primary_to_inertial))) and times ($n) must have the same length"
        ))
        new{R}(times, r_sun, r_secondary, R_primary_to_secondary, R_primary_to_inertial)
    end
end

# Temperature only: R_primary_to_inertial = nothing
BinaryAsteroidEphemerides(times, r_sun, r_secondary, R_primary_to_secondary) =
    BinaryAsteroidEphemerides{Nothing}(times, r_sun, r_secondary, R_primary_to_secondary, nothing)

# With rotation matrices for force/torque computation
BinaryAsteroidEphemerides(times, r_sun, r_secondary, R_primary_to_secondary, R_primary_to_inertial) =
    BinaryAsteroidEphemerides{typeof(R_primary_to_inertial)}(
        times, r_sun, r_secondary, R_primary_to_secondary, R_primary_to_inertial,
    )
