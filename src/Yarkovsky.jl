
# """

# - `a`     : Acceleration
# - `u`     : Eccentric anomaly
# - `orbit` : `OrbitalElements`

# # Return
# - `ȧ` : Drift rate in semi-major axis
# """
# function drift_rate_ȧ(a::AbstractVector, u::Real, orbit::OrbitalElements)
#     @unpack e, n = orbit
#     p = orbit.a * (1 - e^2)  # Semilatus rectum
#     f = u2ν(u, orbit)        # Eccentric anomaly -> True anomaly
#     r = get_r(orbit, u)      # Position @orbital plane frame

#     r = rotateZ(SVector{3}(r), f)  # Rotating frame
#     a = rotateZ(SVector{3}(a), f)  # [a_x, a_y, a_z] -> [a_R, a_T, a_N]

#     ȧ = 2 / n / √(1 - e^2) * (e * sin(f) * a[1] + p / norm(r) * a[2])
# end


