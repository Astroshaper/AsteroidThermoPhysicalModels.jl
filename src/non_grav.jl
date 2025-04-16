

# ****************************************************************
#                  Photon recoil force / torque
# ****************************************************************

"""
    update_thermal_force!(stpm::SingleAsteroidThermoPhysicalModel)

Calculate the thermal force and torque on every face and integrate them over all faces.

# Arguments
- `stpm` : Thermophysical model for a single asteroid
"""
function update_thermal_force!(stpm::SingleAsteroidThermoPhysicalModel)
    stpm.force  .= 0.
    stpm.torque .= 0.

    for i in eachindex(stpm.shape.faces)
        rᵢ = stpm.shape.face_centers[i]
        r̂ᵢ = normalize(rᵢ)
        n̂ᵢ = stpm.shape.face_normals[i]
        aᵢ = stpm.shape.face_areas[i]

        F_sun  = stpm.flux[i, 1]
        F_scat = stpm.flux[i, 2]
        F_rad  = stpm.flux[i, 3]
        Tᵢ = stpm.temperature[begin, i]

        ## Total emittance from face i , Eᵢ [W/m²].
        ## Note that both scattered light and thermal radiation are assumed to be isotropic.
        R_vis = stpm.thermo_params.reflectance_vis[i]
        R_ir  = stpm.thermo_params.reflectance_ir[i]
        ε     = stpm.thermo_params.emissivity[i]
        Eᵢ    = R_vis * F_sun + R_vis * F_scat + R_ir * F_rad + ε * σ_SB * Tᵢ^4

        ## Thermal force on each face
        stpm.face_forces[i] = - 2/3 * Eᵢ * aᵢ / c₀ * n̂ᵢ      # The first term normal to the face
        for visiblefacet in stpm.shape.visiblefacets[i]
            fᵢⱼ = visiblefacet.f
            d̂ᵢⱼ = visiblefacet.d̂

            stpm.face_forces[i] += Eᵢ * aᵢ / c₀ * fᵢⱼ * d̂ᵢⱼ  # The second term of self-heating
        end

        ## Thermal force on the entire shape
        dfᵢ = stpm.face_forces[i]
        stpm.force  .+= (r̂ᵢ ⋅ dfᵢ) * r̂ᵢ  # Thermal force
        stpm.torque .+= rᵢ × dfᵢ         # Thermal torque
    end
end


"""
    update_thermal_force!(btpm::BinaryAsteroidThermoPhysicalModel)

Calculate the thermal force and torque on every face and integrate them over all faces.

# Arguments
- `btpm` : Thermophysical model for a binary asteroid
"""
function update_thermal_force!(btpm::BinaryAsteroidThermoPhysicalModel)
    update_thermal_force!(btpm.pri)
    update_thermal_force!(btpm.sec)
end

