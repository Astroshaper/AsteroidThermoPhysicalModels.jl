#=
grid_params.jl

Numerical grid settings for 1D heat conduction.
=#

"""
    struct GridParams

Numerical grid settings for the 1D heat conduction equation,
shared across all facets.

# Fields
- `z_max`   : Depth of the lower boundary [m]
- `n_depth` : Number of depth nodes
- `Δz`      : Depth step width [m]

# Notes
Currently assumes a uniform depth grid. Variable-spacing support (e.g., finer
near the surface) is planned for a future version.
"""
struct GridParams
    z_max   ::Float64
    n_depth ::Int
    Δz      ::Float64
end

"""
    GridParams(; z_max, n_depth)

Construct `GridParams` from keyword arguments, with `Δz` computed automatically as
`z_max / (n_depth - 1)`, placing nodes uniformly at `0, Δz, 2Δz, …, z_max`.

# Keyword Arguments
- `z_max`   : Depth of the lower boundary [m]
- `n_depth` : Number of depth nodes

# Notes
To specify `Δz` explicitly (advanced use), use the positional constructor
`GridParams(z_max, n_depth, Δz)`.
"""
GridParams(; z_max, n_depth) = GridParams(z_max, n_depth, z_max / (n_depth - 1))
