#=
non_grav.jl

Non-gravitational force calculations for asteroids.
This file implements the thermal recoil effects:
- Yarkovsky effect: Orbital perturbation due to asymmetric thermal emission
- YORP effect: Rotational perturbation due to asymmetric thermal emission
These effects arise from the recoil momentum of photons emitted from the surface.
=#

# ╔═══════════════════════════════════════════════════════════════════╗
# ║                     Photon recoil force / torque                  ║
# ╚═══════════════════════════════════════════════════════════════════╝

"""
    update_thermal_force!(stpm::SingleAsteroidTPM)

Calculate the thermal recoil force (Yarkovsky effect) and torque (YORP effect) on the asteroid
by integrating photon momentum from thermal emission and reflection over all surface facets.

# Arguments
- `stpm::SingleAsteroidTPM` : Thermophysical model for a single asteroid

# Physics
The function calculates non-gravitational effects caused by anisotropic photon emission:
- **Yarkovsky effect**: Net force due to thermal lag causing asymmetric emission
- **YORP effect**: Net torque changing the asteroid's rotation state

# Algorithm
For each facet i, the thermal force is computed as:
```
F_i = -(2/3) × (E_i × A_i)/c × n̂_i + Σⱼ (E_i × A_i)/c × f_ij × d̂_ij
```
where:
- E_i = total emittance from facet i (reflection + thermal emission) [W/m²]
- A_i = area of facet i [m²]
- c = speed of light [m/s]
- n̂_i = outward normal vector of facet i
- f_ij = view factor from facet i to j
- d̂_ij = unit vector from facet i to j

The first term represents direct photon recoil normal to the surface.
The second term accounts for photons intercepted by other facets (self-heating contribution).

# Outputs (stored in stpm)
- `stpm.face_forces` : Thermal force vector on each facet [N]
- `stpm.force` : Total thermal force in body-fixed frame [N]
- `stpm.torque` : Total thermal torque in body-fixed frame [N⋅m]

# Physical Significance
- The force causes orbital drift (Yarkovsky effect)
- The torque changes rotation period and obliquity (YORP effect)
- Both effects are crucial for asteroid orbital evolution

# References
- Bottke Jr, W. F., et al. (2006). The Yarkovsky and YORP effects
- Rozitis, B., & Green, S. F. (2012). The influence of rough surface thermal-infrared beaming
"""
function update_thermal_force!(stpm::SingleAsteroidTPM)
    stpm.force  .= 0.
    stpm.torque .= 0.

    for i in eachindex(stpm.shape.faces)
        rᵢ = stpm.shape.face_centers[i]
        r̂ᵢ = normalize(rᵢ)
        n̂ᵢ = stpm.shape.face_normals[i]
        aᵢ = stpm.shape.face_areas[i]

        F_sun  = stpm.flux_sun[i]
        F_scat = stpm.flux_scat[i]
        F_rad  = stpm.flux_rad[i]
        Tᵢ = stpm.temperature[begin, i]

        ## Total emittance from face i , Eᵢ [W/m²].
        ## Note that both scattered light and thermal radiation are assumed to be isotropic.
        R_vis = stpm.thermo_params.reflectance_vis[i]
        R_ir  = stpm.thermo_params.reflectance_ir[i]
        ε     = stpm.thermo_params.emissivity[i]
        Eᵢ    = R_vis * F_sun + R_vis * F_scat + R_ir * F_rad + ε * σ_SB * Tᵢ^4

        ## Thermal force on each face
        # Photon recoil force: F = -momentum flux = -Energy flux / c
        # The factor 2/3 comes from Lambertian emission (isotropic in hemisphere)
        # For Lambertian surface: ∫cos(θ)dΩ = 2π/3 over hemisphere
        stpm.face_forces[i] = - 2/3 * Eᵢ * aᵢ / c₀ * n̂ᵢ  # Direct recoil force normal to face
        for visiblefacet in stpm.shape.visiblefacets[i]
            fᵢⱼ = visiblefacet.f
            d̂ᵢⱼ = visiblefacet.d̂

            # Self-heating contribution: photons re-absorbed by other faces
            # No 2/3 factor here because the direction is already specified by d̂ᵢⱼ
            stpm.face_forces[i] += Eᵢ * aᵢ / c₀ * fᵢⱼ * d̂ᵢⱼ  # Re-absorption recoil force
        end

        ## Thermal force on the entire shape
        dfᵢ = stpm.face_forces[i]
        stpm.force  .+= (r̂ᵢ ⋅ dfᵢ) * r̂ᵢ  # Thermal force
        stpm.torque .+= rᵢ × dfᵢ         # Thermal torque
    end
end


"""
    update_thermal_force!(btpm::BinaryAsteroidTPM)

Calculate the thermal force and torque on every face and integrate them over all faces.

# Arguments
- `btpm` : Thermophysical model for a binary asteroid
"""
function update_thermal_force!(btpm::BinaryAsteroidTPM)
    update_thermal_force!(btpm.pri)
    update_thermal_force!(btpm.sec)
end
