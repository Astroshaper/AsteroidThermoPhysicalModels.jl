

# These functions will be used if the surface roughness models are implemented
# into the thermophysical simulation.

################################################################
#                  Concave spherical segment
################################################################

# crater_curvature_radius(r, h) = (r^2 + h^2) / 2h

# function concave_spherical_segment(r, h, xc, yc, x, y)
#     d² = (x - xc)^2 + (y - yc)^2
#     d = √d²  # Distance from the crater center

#     if d > r
#         z = 0.
#     else
#         R = crater_curvature_radius(r, h)
#         z = R - h - √(R^2 - d²)
#     end
#     z
# end

# function concave_spherical_segment(r, h; Nx=2^5, Ny=2^5, xc=0.5, yc=0.5)
#     xs = LinRange(0, 1, Nx + 1)
#     ys = LinRange(0, 1, Ny + 1)
#     zs = [concave_spherical_segment(r, h, xc, yc, x, y) for x in xs, y in ys]

#     xs, ys, zs
# end

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