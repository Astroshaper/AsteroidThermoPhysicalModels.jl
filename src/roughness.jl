

################################################################
#                  Concave spherical segment
################################################################

"""
    crater_curvature_radius(r, h) -> R

Calculate the curvature radius of a concave spherical segment (crater shape).

# Arguments
- `r::Real`: Crater radius [m]
- `h::Real`: Crater depth [m]

# Returns
- `R::Real`: Curvature radius of the spherical cap [m]

# Mathematical Formula
The curvature radius `R` of a spherical cap is given by:
```math
R = \\frac{r^2 + h^2}{2h}
```
where `r` is the crater rim radius and `h` is the crater depth.

# Example
```julia
# A crater with 10m radius and 2m depth
R = crater_curvature_radius(10.0, 2.0)  # Returns 26.0
```

# Notes
- The crater is modeled as a spherical cap (portion of a sphere)
- Larger R values correspond to shallower craters
- As h → 0, R → ∞ (flat surface)
"""
crater_curvature_radius(r::Real, h::Real) = (r^2 + h^2) / 2h

"""
    concave_spherical_segment(r, h, xc, yc, x, y) -> z

Calculate the z-coordinate (depth) of a point on a concave spherical segment (crater).

# Arguments
- `r::Real` : Crater radius [m]
- `h::Real` : Crater depth [m]
- `xc::Real`: x-coordinate of crater center [m]
- `yc::Real`: y-coordinate of crater center [m]
- `x::Real` : x-coordinate where to calculate z [m]
- `y::Real` : y-coordinate where to calculate z [m]

# Returns
- `z::Real`: Depth below the reference surface (negative value) [m]

# Mathematical Model
The crater is modeled as a spherical cap. For a point (x,y) within the crater:
```math
z = R - h - \\sqrt{R^2 - d^2}
```
where:
- `R` is the curvature radius calculated from `crater_curvature_radius(r, h)`
- `d` is the distance from the crater center: `d = √((x-xc)² + (y-yc)²)`
- `z = 0` outside the crater (when `d > r`)

# Example
```julia
# A crater at (0.5, 0.5) with radius 0.2 and depth 0.05
z = concave_spherical_segment(0.2, 0.05, 0.5, 0.5, 0.6, 0.5)  # Point at crater rim
```

# Notes
- Returns `z ≤ 0` (crater is a depression)
- Returns `z = 0` for points outside the crater
- The reference surface is at `z = 0`
"""
function concave_spherical_segment(r::Real, h::Real, xc::Real, yc::Real, x::Real, y::Real)
    d² = (x - xc)^2 + (y - yc)^2
    d = √d²  # Distance from the crater center

    if d > r
        z = 0.
    else
        R = crater_curvature_radius(r, h)
        z = R - h - √(R^2 - d²)
    end
    z
end

"""
    concave_spherical_segment(r, h; Nx=2^5, Ny=2^5, xc=0.5, yc=0.5) -> xs, ys, zs

Generate a grid representation of a concave spherical segment (crater) on a unit square domain.

# Arguments
- `r::Real` : Crater radius (in normalized coordinates, 0-1)
- `h::Real` : Crater depth (in same units as radius)

# Keyword Arguments  
- `Nx::Integer=32`: Number of grid points in the x-direction
- `Ny::Integer=32`: Number of grid points in the y-direction
- `xc::Real=0.5`: x-coordinate of crater center (normalized, 0-1)
- `yc::Real=0.5`: y-coordinate of crater center (normalized, 0-1)

# Returns
- `xs::LinRange`: x-coordinates of grid points (0 to 1)
- `ys::LinRange`: y-coordinates of grid points (0 to 1)
- `zs::Matrix`: z-coordinates (depths) at each grid point

# Example
```julia
# Create a crater grid with radius 0.3 and depth 0.1
xs, ys, zs = concave_spherical_segment(0.3, 0.1; Nx=64, Ny=64)

# Off-center crater
xs, ys, zs = concave_spherical_segment(0.2, 0.05; xc=0.3, yc=0.7)
```

# Notes
- The domain is normalized to [0,1] × [0,1]
- The crater parameters should be scaled accordingly
- Grid has (Nx+1) × (Ny+1) total points
- Useful for generating synthetic rough surfaces
"""
function concave_spherical_segment(r::Real, h::Real; Nx::Integer=2^5, Ny::Integer=2^5, xc::Real=0.5, yc::Real=0.5)
    xs = LinRange(0, 1, Nx + 1)
    ys = LinRange(0, 1, Ny + 1)
    zs = [concave_spherical_segment(r, h, xc, yc, x, y) for x in xs, y in ys]

    xs, ys, zs
end


################################################################
#                 Parallel sinusoidal trenches
################################################################

# function parallel_sinusoidal_trenches(Ct, N_trench, x)
#     z = Ct + Ct * sin(2π * (N_trench + 0.5) * x)
# end

# function parallel_sinusoidal_trenches(Ct, N_trench; Nx=2^5, Ny=2^5)
#     xs = LinRange(0, 1, Nx + 1)
#     ys = LinRange(0, 1, Ny + 1)
#     zs = [parallel_sinusoidal_trenches(Ct, N_trench, x) for x in xs, y in ys]

#     xs, ys, zs 
# end


################################################################
#                   Random Gaussian surface
################################################################

# TODO: Functions to generate random Gaussian surface will be implmented.


################################################################
#                       Fractal surface
################################################################

# TODO: Functions to generate fractal surface will be implmented.

