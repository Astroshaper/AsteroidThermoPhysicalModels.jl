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
    update_thermal_force!(state::SingleAsteroidThermoPhysicalState)

Calculate the thermal recoil force (Yarkovsky effect) and torque (YORP effect) on the asteroid
by integrating photon momentum from thermal emission and reflection over all surface facets.

# Arguments
- `state::SingleAsteroidThermoPhysicalState` : Thermophysical simulation state for a single asteroid

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
The second term accounts for photons intercepted by other facets: their momentum stays with
the body, so it cancels part of the recoil. It is the momentum counterpart of the energy
re-absorbed in self-heating, and is therefore applied only when `with_self_heating` is
enabled; without it, every emitted photon counts as having left the body.

# Outputs (stored in state)
- `state.face_forces` : Thermal force vector on each facet [N]
- `state.force` : Net thermal force `Σᵢ F_i` in the body-fixed frame [N]
- `state.torque` : Net thermal torque `Σᵢ r_i × F_i` about the body-fixed origin [N⋅m]; the
  origin is assumed to be the centre of mass

# Physical Significance
- The force causes orbital drift (Yarkovsky effect)
- The torque changes rotation period and obliquity (YORP effect)
- Both effects are crucial for asteroid orbital evolution

# References
- Bottke Jr, W. F., et al. (2006). The Yarkovsky and YORP effects
- Rozitis, B., & Green, S. F. (2012). The influence of rough surface thermal-infrared beaming
"""
function update_thermal_force!(state::SingleAsteroidThermoPhysicalState)
    state.force  .= 0.
    state.torque .= 0.

    with_self_heating = state.problem.with_self_heating
    with_self_heating && _require_face_visibility_graph(state.problem.shape, "with_self_heating")

    for i in eachindex(state.problem.shape.faces)
        rᵢ = state.problem.shape.face_centers[i]
        n̂ᵢ = state.problem.shape.face_normals[i]
        aᵢ = state.problem.shape.face_areas[i]

        F_sun  = state.flux_sun[i]
        F_scat = state.flux_scat[i]
        F_rad  = state.flux_rad[i]
        Tᵢ = state.temperature[begin, i]

        ## Total emittance from face i , Eᵢ [W/m²].
        ## Note that both scattered light and thermal radiation are assumed to be isotropic.
        R_vis = state.problem.thermo_params.reflectance_vis[i]
        R_ir  = state.problem.thermo_params.reflectance_ir[i]
        ε     = state.problem.thermo_params.emissivity[i]
        Eᵢ    = R_vis * F_sun + R_vis * F_scat + R_ir * F_rad + ε * σ_SB * Tᵢ^4

        ## Thermal force on each face
        # Photon recoil force: F = -momentum flux = -Energy flux / c
        # The factor 2/3 comes from Lambertian emission (isotropic in hemisphere)
        # For Lambertian surface: ∫cos(θ)dΩ = 2π/3 over hemisphere
        state.face_forces[i] = - 2/3 * Eᵢ * aᵢ / c₀ * n̂ᵢ  # Direct recoil force normal to face
        
        if with_self_heating
            # Face properties visible from `i`: View factors and directions
            view_factors = get_view_factors(state.problem.shape.face_visibility_graph, i)
            directions = get_visible_face_directions(state.problem.shape.face_visibility_graph, i)
            
            for (fᵢⱼ, d̂ᵢⱼ) in zip(view_factors, directions)
                # Self-heating contribution: photons re-absorbed by other faces
                # No 2/3 factor here because the direction is already specified by d̂ᵢⱼ
                state.face_forces[i] += Eᵢ * aᵢ / c₀ * fᵢⱼ * d̂ᵢⱼ  # Re-absorption recoil force
            end
        end

        ## Thermal force and torque on the entire shape. The net force is the plain sum of
        ## the face forces: where a force acts does not enter the motion of the centre of
        ## mass, only the torque. The torque is taken about the body-fixed origin, which is
        ## assumed to be the centre of mass.
        dfᵢ = state.face_forces[i]
        state.force  .+= dfᵢ        # Thermal force
        state.torque .+= rᵢ × dfᵢ   # Thermal torque
    end
end


"""
    update_thermal_force!(state::BinaryAsteroidThermoPhysicalState)

Calculate the thermal force and torque on every face and integrate them over all faces.

# Arguments
- `state` : Thermophysical simulation state for a binary asteroid
"""
function update_thermal_force!(state::BinaryAsteroidThermoPhysicalState)
    update_thermal_force!(state.primary)
    update_thermal_force!(state.secondary)
end
