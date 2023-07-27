

# ****************************************************************
#                  Photon recoil force / torque
# ****************************************************************

"""
    update_thermal_force!(shape::ShapeModel, params)
    update_thermal_force!(shape::ShapeModel, A_B, ε)

Calculate the thermal force and torque on every facet and update their totals.
"""
update_thermal_force!(shape::ShapeModel, params::AbstractThermoParams) = update_thermal_force!(shape, params.A_B, params.ε)

function update_thermal_force!(shape::ShapeModel, A_B, ε)
    shape.force  .= 0
    shape.torque .= 0
    for i in eachindex(shape.faces)
        rᵢ = shape.face_centers[i]
        n̂ᵢ = shape.face_normals[i]
        aᵢ = shape.face_areas[i]

        F_scat = shape.flux[i, 2]
        Tᵢ = shape.facets[i].temps[begin]

        ## Total amount of scattered light and radiation [W/m²].
        ## Note that both are assumed to be isotropic.
        Eᵢ = (A_B isa Real ? A_B : A_B[i]) * F_scat + (ε isa Real ? ε : ε[i]) * σ_SB * Tᵢ^4

        ## Thermal force on each face
        @. shape.facets[i].force = - 2/3 * Eᵢ * aᵢ / c₀ * n̂ᵢ      # The first term normal to the face
        for visiblefacet in shape.facets[i].visiblefacets
            fᵢⱼ = visiblefacet.f
            d̂ᵢⱼ = visiblefacet.d̂

            @. shape.facets[i].force += Eᵢ * aᵢ / c₀ * fᵢⱼ * d̂ᵢⱼ  # The second term of self-heating
        end

        ## Thermal force on the entire shape
        r̂ᵢ = normalize(rᵢ)
        dfᵢ = SVector{3}(shape.facets[i].force)

        shape.force  .+= (r̂ᵢ ⋅ dfᵢ) * r̂ᵢ  # Thermal force
        shape.torque .+= rᵢ × dfᵢ         # Thermal torque
    end
end

