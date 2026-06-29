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
- `Δz`      : Depth step width [m]
- `n_depth` : Number of depth nodes

# Notes
Currently assumes a uniform depth grid. Variable-spacing support (e.g., finer
near the surface) is planned for a future version.
"""
struct GridParams
    z_max   ::Float64
    Δz      ::Float64
    n_depth ::Int
end
