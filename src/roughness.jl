

################################################################
#                  Concave spherical segment
################################################################

"""
    crater_curvature_radius(r, h) -> R

Return the curvature radius of a concave spherical segment.

# Arguments
- `r`: Crater radius
- `h`: Crater depth
"""
crater_curvature_radius(r::Real, h::Real) = (r^2 + h^2) / 2h

"""
    concave_spherical_segment(r, h, xc, yc, x, y) -> z

Return the z-coordinate of a concave spherical segment.

# Arguments
- `r` : Crater radius
- `h` : Crater depth
- `xc`: x-coordinate of crater center
- `yc`: y-coordinate of crater center
- `x` : x-coordinate where to calculate z
- `y` : y-coordinate where to calculate z
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

Return (x, y, z) grid of a concave spherical segment.

# Arguments
- `r` : Crater radius
- `h` : Crater depth
- `Nx`: Number of nodes in the x-direction
- `Ny`: Number of nodes in the y-direction
- `xc`: x-coordinate of crater center
- `yc`: y-coordinate of crater center
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

# TO DO: Functions to generate random Gaussian surface will be implmented.


################################################################
#                       Fractal surface
################################################################

# TO DO: Functions to generate fractal surface will be implmented.

