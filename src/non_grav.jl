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

# Total emittance of face `i`, Eᵢ [W/m²]: reflected sunlight and scattered light, reflected
# thermal radiation, and thermal emission. Both reflection and emission are taken as isotropic.
function _emittance(state::SingleAsteroidThermoPhysicalState, i::Integer)
    R_vis = state.problem.thermo_params.reflectance_vis[i]
    R_ir  = state.problem.thermo_params.reflectance_ir[i]
    ε     = state.problem.thermo_params.emissivity[i]
    Tᵢ    = state.temperature[begin, i]

    R_vis * state.flux_sun[i] + R_vis * state.flux_scat[i] + R_ir * state.flux_rad[i] + ε * σ_SB * Tᵢ^4
end


# Recoil force on every face of `shape`, written to `state.face_forces`. The net force and
# torque are left to `_accumulate_force_torque!`, so that a caller may adjust individual face
# forces in between. The shape is passed explicitly, as for the flux updates; the function
# is reused for the sub-face states of the roughness models.
function _update_face_forces!(state::SingleAsteroidThermoPhysicalState, shape::ShapeModel)
    with_self_heating = state.problem.with_self_heating
    with_self_heating && _require_face_visibility_graph(shape, "with_self_heating")

    for i in eachindex(shape.faces)
        n̂ᵢ = shape.face_normals[i]
        aᵢ = shape.face_areas[i]
        Eᵢ = _emittance(state, i)

        ## Thermal force on each face
        # Photon recoil force: F = -momentum flux = -Energy flux / c
        # The factor 2/3 comes from Lambertian emission (isotropic in hemisphere)
        # For Lambertian surface: ∫cos(θ)dΩ = 2π/3 over hemisphere
        state.face_forces[i] = - 2/3 * Eᵢ * aᵢ / c₀ * n̂ᵢ  # Direct recoil force normal to face

        if with_self_heating
            # Face properties visible from `i`: View factors and directions
            view_factors = get_view_factors(shape.face_visibility_graph, i)
            directions = get_visible_face_directions(shape.face_visibility_graph, i)

            for (fᵢⱼ, d̂ᵢⱼ) in zip(view_factors, directions)
                # Self-heating contribution: photons re-absorbed by other faces
                # No 2/3 factor here because the direction is already specified by d̂ᵢⱼ
                state.face_forces[i] += Eᵢ * aᵢ / c₀ * fᵢⱼ * d̂ᵢⱼ  # Re-absorption recoil force
            end
        end
    end
end


# Net force and torque from `state.face_forces`. The net force is the plain sum of the face
# forces: where a force acts does not enter the motion of the centre of mass, only the torque.
# The torque is taken about the body-fixed origin, which is assumed to be the centre of mass.
function _accumulate_force_torque!(state::SingleAsteroidThermoPhysicalState, shape::ShapeModel)
    state.force  .= 0.
    state.torque .= 0.

    for i in eachindex(shape.faces)
        rᵢ  = shape.face_centers[i]
        dfᵢ = state.face_forces[i]
        state.force  .+= dfᵢ        # Thermal force
        state.torque .+= rᵢ × dfᵢ   # Thermal torque
    end
end


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

# Surface roughness
When the shape carries surface roughness, a face with a roughness model gets its recoil from
the sub-faces of its model instead of the smooth-surface baseline: the sub-face forces are
summed in the local frame, scaled from the patch to the parent's area by `Aᵢ / A_proj`, and
rotated into the body frame. With `with_self_heating`, the photons that escape the patch
towards the sky and are intercepted by other faces are also accounted for (isotropic patch
emission). The torque takes the centre of the parent face as the point of action; the torque
of the patch about its own centre is neglected. The sub-face states' own `force` and `torque`
are not used. Only the derived quantities of the faces are overwritten; their fluxes and
temperatures remain the smooth-surface solution.

# References
- Bottke Jr, W. F., et al. (2006). The Yarkovsky and YORP effects
- Rozitis, B., & Green, S. F. (2012). The influence of rough surface thermal-infrared beaming
"""
function update_thermal_force!(state::SingleAsteroidThermoPhysicalState)
    shape             = state.problem.shape
    with_self_heating = state.problem.with_self_heating

    _update_face_forces!(state, shape)  # smooth-surface baseline on every face

    # Sub-face level, for every face `i` that carries a roughness model (empty loop for a
    # smooth surface): the recoil on each sub-face `j` in the local frame of the model,
    # including the momentum of the photons intercepted by the other sub-faces (self-heating
    # inside a roughness model is always on). The roughness model is a patch that represents
    # the surface of the parent face statistically, so the sum is counted for the parent's
    # area rather than for the area of the patch, and rotated into the body frame:
    #
    #     F_i = (A_i / A_proj) × R_iᵀ Σⱼ f_j,    A_proj = Σⱼ a_j (n̂_j ⋅ ẑ)⁺
    #
    # where `f_j` and `a_j` are the force and area of sub-face `j` in the units of the model,
    # `A_proj` is the area of the model projected onto its reference plane (`projected_area`
    # of AsteroidShapeModels.jl), and `R_i` rotates from the body frame to the local frame.
    # The `scale` of the roughness model cancels: the force on one patch grows as `scale²`
    # and the number of patches covering the face falls as `scale⁻²`. The result replaces
    # `face_forces[i]`.
    for (i, k) in enumerate(state.face_roughness_indices)
        k == 0 && continue
        rs = state.roughness_states[k]
        roughness_shape = rs.problem.shape
        update_thermal_force!(rs)  # includes the re-absorption inside the roughness model

        # Representative patch: count the model-unit forces for the parent's area. The
        # rotation is linear, so the sub-face forces are summed first and rotated once.
        patches_per_face = shape.face_areas[i] / projected_area(roughness_shape)
        f_sum = sum(rs.face_forces)
        state.face_forces[i] = patches_per_face * transform_physical_vector_local_to_global(shape, i, f_sum)

        # Photons that leave the roughness model towards the sky and are intercepted by other
        # faces: the momentum counterpart of the external irradiation, isotropic. The power
        # escaping the model, P_sky = (A_i / A_proj) Σⱼ E_j a_j f_sky,j, contributes
        # P_sky / c × Σₖ f_ik d̂_ik.
        if with_self_heating
            graph_sub = roughness_shape.face_visibility_graph
            P_sky = patches_per_face * sum(eachindex(roughness_shape.faces)) do j
                _emittance(rs, j) * roughness_shape.face_areas[j] * _sky_view_factor(graph_sub, j)
            end
            view_factors = get_view_factors(shape.face_visibility_graph, i)
            directions   = get_visible_face_directions(shape.face_visibility_graph, i)
            for (fᵢₖ, d̂ᵢₖ) in zip(view_factors, directions)
                state.face_forces[i] += P_sky / c₀ * fᵢₖ * d̂ᵢₖ
            end
        end
    end

    _accumulate_force_torque!(state, shape)
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
