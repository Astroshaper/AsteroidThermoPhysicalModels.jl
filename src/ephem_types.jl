#=
ephem_types.jl

Formal ephemerides types for thermophysical simulations.
These replace the informal NamedTuple previously used as `ephem`.

Type hierarchy:
    AbstractAsteroidEphemerides
    ├── AbstractSingleAsteroidEphemerides
    │   ├── SingleAsteroidEphemerides
    │   └── SingleAsteroidEphemeridesWithDynamics
    └── AbstractBinaryAsteroidEphemerides
        ├── BinaryAsteroidEphemerides
        └── BinaryAsteroidEphemeridesWithDynamics
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
    struct SingleAsteroidEphemerides <: AbstractSingleAsteroidEphemerides

Ephemerides for a single-asteroid thermophysical simulation.
Only temperature is computed; force and torque are not.

# Fields
- `times` : Simulation timesteps [s]
- `r_sun` : Sun position vector in the body-fixed frame [m]
"""
struct SingleAsteroidEphemerides <: AbstractSingleAsteroidEphemerides
    times ::Vector{Float64}
    r_sun ::Vector{SVector{3, Float64}}

    function SingleAsteroidEphemerides(times, r_sun)
        length(times) == length(r_sun) || throw(DimensionMismatch(
            "times ($(length(times))) and r_sun ($(length(r_sun))) must have the same length"
        ))
        new(times, r_sun)
    end
end

"""
    struct SingleAsteroidEphemeridesWithDynamics <: AbstractSingleAsteroidEphemerides

Ephemerides for a single-asteroid thermophysical simulation with force and torque output
in the inertial frame. Requires the body orientation at each timestep.

# Fields
- `times`              : Simulation timesteps [s]
- `r_sun`              : Sun position vector in the body-fixed frame [m]
- `R_body_to_inertial` : Passive rotation matrix from the body-fixed frame to the inertial frame
"""
struct SingleAsteroidEphemeridesWithDynamics <: AbstractSingleAsteroidEphemerides
    times              ::Vector{Float64}
    r_sun              ::Vector{SVector{3, Float64}}
    R_body_to_inertial ::Vector{SMatrix{3,3,Float64,9}}

    function SingleAsteroidEphemeridesWithDynamics(times, r_sun, R_body_to_inertial)
        n = length(times)
        length(r_sun)              == n || throw(DimensionMismatch(
            "r_sun ($(length(r_sun))) and times ($n) must have the same length"
        ))
        length(R_body_to_inertial) == n || throw(DimensionMismatch(
            "R_body_to_inertial ($(length(R_body_to_inertial))) and times ($n) must have the same length"
        ))
        new(times, r_sun, R_body_to_inertial)
    end
end

"""
    SingleAsteroidEphemeridesWithDynamics(ephem, R_body_to_inertial)

Upgrade constructor: extend a `SingleAsteroidEphemerides` with body orientation data.
"""
function SingleAsteroidEphemeridesWithDynamics(
    ephem              ::SingleAsteroidEphemerides,
    R_body_to_inertial ::Vector{SMatrix{3,3,Float64,9}},
)
    SingleAsteroidEphemeridesWithDynamics(ephem.times, ephem.r_sun, R_body_to_inertial)
end

# ╔═══════════════════════════════════════════════════════════════════╗
# ║                     Binary asteroid                               ║
# ╚═══════════════════════════════════════════════════════════════════╝

"""
    struct BinaryAsteroidEphemerides <: AbstractBinaryAsteroidEphemerides

Ephemerides for a binary-asteroid thermophysical simulation.
Only temperature is computed; force and torque are not.

# Fields
- `times`                  : Simulation timesteps [s]
- `r_sun`                  : Sun position vector in the primary body-fixed frame [m]
- `r_secondary`            : Secondary position vector in the primary body-fixed frame [m]
- `R_primary_to_secondary` : Passive rotation matrix from the primary body-fixed frame to the secondary body-fixed frame
"""
struct BinaryAsteroidEphemerides <: AbstractBinaryAsteroidEphemerides
    times                  ::Vector{Float64}
    r_sun                  ::Vector{SVector{3, Float64}}
    r_secondary            ::Vector{SVector{3, Float64}}
    R_primary_to_secondary ::Vector{SMatrix{3,3,Float64,9}}

    function BinaryAsteroidEphemerides(times, r_sun, r_secondary, R_primary_to_secondary)
        n = length(times)
        length(r_sun)                  == n || throw(DimensionMismatch(
            "r_sun ($(length(r_sun))) and times ($n) must have the same length"
        ))
        length(r_secondary)            == n || throw(DimensionMismatch(
            "r_secondary ($(length(r_secondary))) and times ($n) must have the same length"
        ))
        length(R_primary_to_secondary) == n || throw(DimensionMismatch(
            "R_primary_to_secondary ($(length(R_primary_to_secondary))) and times ($n) must have the same length"
        ))
        new(times, r_sun, r_secondary, R_primary_to_secondary)
    end
end

"""
    struct BinaryAsteroidEphemeridesWithDynamics <: AbstractBinaryAsteroidEphemerides

Ephemerides for a binary-asteroid thermophysical simulation with force and torque output
in the inertial frame. Requires the primary body orientation at each timestep.

The secondary-to-inertial rotation can be derived as
`R_primary_to_inertial[i] * R_primary_to_secondary[i]'` when needed.

# Fields
- `times`                  : Simulation timesteps [s]
- `r_sun`                  : Sun position vector in the primary body-fixed frame [m]
- `r_secondary`            : Secondary position vector in the primary body-fixed frame [m]
- `R_primary_to_secondary` : Passive rotation matrix from the primary body-fixed frame to the secondary body-fixed frame
- `R_primary_to_inertial`  : Passive rotation matrix from the primary body-fixed frame to the inertial frame
"""
struct BinaryAsteroidEphemeridesWithDynamics <: AbstractBinaryAsteroidEphemerides
    times                  ::Vector{Float64}
    r_sun                  ::Vector{SVector{3, Float64}}
    r_secondary            ::Vector{SVector{3, Float64}}
    R_primary_to_secondary ::Vector{SMatrix{3,3,Float64,9}}
    R_primary_to_inertial  ::Vector{SMatrix{3,3,Float64,9}}

    function BinaryAsteroidEphemeridesWithDynamics(
        times, r_sun, r_secondary, R_primary_to_secondary, R_primary_to_inertial,
    )
        n = length(times)
        length(r_sun)                  == n || throw(DimensionMismatch(
            "r_sun ($(length(r_sun))) and times ($n) must have the same length"
        ))
        length(r_secondary)            == n || throw(DimensionMismatch(
            "r_secondary ($(length(r_secondary))) and times ($n) must have the same length"
        ))
        length(R_primary_to_secondary) == n || throw(DimensionMismatch(
            "R_primary_to_secondary ($(length(R_primary_to_secondary))) and times ($n) must have the same length"
        ))
        length(R_primary_to_inertial)  == n || throw(DimensionMismatch(
            "R_primary_to_inertial ($(length(R_primary_to_inertial))) and times ($n) must have the same length"
        ))
        new(times, r_sun, r_secondary, R_primary_to_secondary, R_primary_to_inertial)
    end
end

"""
    BinaryAsteroidEphemeridesWithDynamics(ephem, R_primary_to_inertial)

Upgrade constructor: extend a `BinaryAsteroidEphemerides` with primary body orientation data.
"""
function BinaryAsteroidEphemeridesWithDynamics(
    ephem                 ::BinaryAsteroidEphemerides,
    R_primary_to_inertial ::Vector{SMatrix{3,3,Float64,9}},
)
    BinaryAsteroidEphemeridesWithDynamics(
        ephem.times, ephem.r_sun, ephem.r_secondary, ephem.R_primary_to_secondary,
        R_primary_to_inertial,
    )
end
