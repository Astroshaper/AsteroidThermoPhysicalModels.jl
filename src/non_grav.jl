

# ****************************************************************
#                  Photon recoil force / torque
# ****************************************************************

"""
    update_thermal_force!(stpm::SingleTPM, nₜ::Integer)

Calculate the thermal force and torque on every face and integrate them over all faces.

# Arguments
- `stpm` : Thermophysical model for a single asteroid
- `nₜ`   : Index of time step
"""
function update_thermal_force!(stpm::SingleTPM, nₜ::Integer)
    stpm.force  .= 0.
    stpm.torque .= 0.

    for i in eachindex(stpm.shape.faces)
        rᵢ = stpm.shape.face_centers[i]
        r̂ᵢ = normalize(rᵢ)
        n̂ᵢ = stpm.shape.face_normals[i]
        aᵢ = stpm.shape.face_areas[i]

        F_scat = stpm.flux[i, 2]
        Tᵢ = stpm.temperature[begin, i, nₜ]

        ## Total amount of scattered light and radiation, Eᵢ [W/m²].
        ## Note that both are assumed to be isotropic.
        A_B = (stpm.thermo_params.A_B isa Real ? stpm.thermo_params.A_B : stpm.thermo_params.A_B[i])
        ε   = (stpm.thermo_params.ε   isa Real ? stpm.thermo_params.ε   : stpm.thermo_params.ε[i])
        Eᵢ  = A_B * F_scat + ε * σ_SB * Tᵢ^4

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
    update_thermal_force!(btpm::BinaryTPM, nₜ::Integer)

Calculate the thermal force and torque on every face and integrate them over all faces.

# Arguments
- `btpm` : Thermophysical model for a binary asteroid
- `nₜ`   : Index of time step
"""
function update_thermal_force!(btpm::BinaryTPM, nₜ::Integer)
    update_thermal_force!(btpm.pri, nₜ)
    update_thermal_force!(btpm.sec, nₜ)
end

