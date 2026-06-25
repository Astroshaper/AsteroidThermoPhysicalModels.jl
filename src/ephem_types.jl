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
# ║                     Conversion helpers                            ║
# ╚═══════════════════════════════════════════════════════════════════╝

# Two methods per helper: the first is a no-op fast path when the element type
# already matches (avoids element-wise allocation); the second converts each
# element so users need not import StaticArrays themselves.

_to_svec3_vec(v::AbstractVector{SVector{3,Float64}}) = convert(Vector{SVector{3,Float64}}, v)
_to_svec3_vec(v::AbstractVector) = [SVector{3,Float64}(x) for x in v]

_to_smat33_vec(v::AbstractVector{SMatrix{3,3,Float64,9}}) = convert(Vector{SMatrix{3,3,Float64,9}}, v)
_to_smat33_vec(v::AbstractVector) = [SMatrix{3,3,Float64,9}(x) for x in v]

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
    SingleAsteroidEphemerides(times, r_sun, nothing)
        -> SingleAsteroidEphemerides{Nothing}
    SingleAsteroidEphemerides(times, r_sun, R_body_to_inertial)
        -> SingleAsteroidEphemerides{Vector{SMatrix{3,3,Float64,9}}}

Plain `AbstractVector` / `AbstractMatrix` elements in `r_sun` and
`R_body_to_inertial` are automatically converted to the corresponding
`SVector{3,Float64}` / `SMatrix{3,3,Float64,9}` types, so importing
`StaticArrays` in user code is not required.
"""
struct SingleAsteroidEphemerides{R <: Union{Nothing, AbstractVector}} <: AbstractSingleAsteroidEphemerides
    times              ::Vector{Float64}
    r_sun              ::Vector{SVector{3, Float64}}
    R_body_to_inertial ::R

    function SingleAsteroidEphemerides{R}(
        times              ::AbstractVector,
        r_sun              ::AbstractVector,
        R_body_to_inertial ::R,
    ) where {R <: Union{Nothing, AbstractVector}}
        n = length(times)
        length(r_sun) == n || throw(DimensionMismatch(
            "r_sun ($(length(r_sun))) and times ($n) must have the same length"
        ))
        R_body_to_inertial === nothing || length(R_body_to_inertial) == n || throw(DimensionMismatch(
            "R_body_to_inertial ($(length(R_body_to_inertial))) and times ($n) must have the same length"
        ))
        new{R}(times, _to_svec3_vec(r_sun), R_body_to_inertial)
    end
end

# Temperature only (R_body_to_inertial = nothing), implicit or explicit.
SingleAsteroidEphemerides(times, r_sun) =
    SingleAsteroidEphemerides{Nothing}(times, r_sun, nothing)
SingleAsteroidEphemerides(times, r_sun, ::Nothing) =
    SingleAsteroidEphemerides{Nothing}(times, r_sun, nothing)

# With rotation matrices for force/torque computation.
# Converted here (before the inner constructor) so that the type parameter R is
# correctly inferred as Vector{SMatrix{3,3,Float64,9}}.
function SingleAsteroidEphemerides(times, r_sun, R_body_to_inertial)
    R_conv = _to_smat33_vec(R_body_to_inertial)
    SingleAsteroidEphemerides{typeof(R_conv)}(times, r_sun, R_conv)
end


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
    BinaryAsteroidEphemerides(times, r_sun, r_secondary, R_primary_to_secondary, nothing)
        -> BinaryAsteroidEphemerides{Nothing}
    BinaryAsteroidEphemerides(times, r_sun, r_secondary, R_primary_to_secondary, R_primary_to_inertial)
        -> BinaryAsteroidEphemerides{Vector{SMatrix{3,3,Float64,9}}}

Plain `AbstractVector` / `AbstractMatrix` elements are automatically converted
to `SVector{3,Float64}` / `SMatrix{3,3,Float64,9}` types in all fields.
"""
struct BinaryAsteroidEphemerides{R <: Union{Nothing, AbstractVector}} <: AbstractBinaryAsteroidEphemerides
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
    ) where {R <: Union{Nothing, AbstractVector}}
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
        new{R}(
            times,
            _to_svec3_vec(r_sun),
            _to_svec3_vec(r_secondary),
            _to_smat33_vec(R_primary_to_secondary),
            R_primary_to_inertial,
        )
    end
end

# Temperature only (R_primary_to_inertial = nothing), implicit or explicit.
BinaryAsteroidEphemerides(times, r_sun, r_secondary, R_primary_to_secondary) =
    BinaryAsteroidEphemerides{Nothing}(times, r_sun, r_secondary, R_primary_to_secondary, nothing)
BinaryAsteroidEphemerides(times, r_sun, r_secondary, R_primary_to_secondary, ::Nothing) =
    BinaryAsteroidEphemerides{Nothing}(times, r_sun, r_secondary, R_primary_to_secondary, nothing)

# With rotation matrices for force/torque computation.
# Converted here (before the inner constructor) so that the type parameter R is
# correctly inferred as Vector{SMatrix{3,3,Float64,9}}.
function BinaryAsteroidEphemerides(times, r_sun, r_secondary, R_primary_to_secondary, R_primary_to_inertial)
    R_conv = _to_smat33_vec(R_primary_to_inertial)
    BinaryAsteroidEphemerides{typeof(R_conv)}(times, r_sun, r_secondary, R_primary_to_secondary, R_conv)
end
