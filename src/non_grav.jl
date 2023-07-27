

# ****************************************************************
#                  Photon recoil force / torque
# ****************************************************************

"""
    update_thermal_force!(shape::ShapeModel, params)
    update_thermal_force!(shape::ShapeModel, A_B, ε)
    update_thermal_force!(facet::Facet, A_B::Real, ε::Real, n̂, a)

Calculate the thermal force and torque on every facet and update the total force/torque.
"""
update_thermal_force!(shape::ShapeModel, params) = update_thermal_force!(shape, params.A_B, params.ε)

function update_thermal_force!(shape::ShapeModel, A_B, ε)
    shape.force  .= 0
    shape.torque .= 0
    for i in eachindex(shape.faces)
        update_thermal_force!(shape.facets[i], (A_B isa Real ? A_B : A_B[i]), (ε isa Real ? ε : ε[i]), shape.face_normals[i], shape.face_areas[i])

        r = shape.face_centers[i]
        r̂ = normalize(r)
        df = SVector{3}(shape.facets[i].force)

        shape.force  .+= (r̂ ⋅ df) * r̂  # Photon recoil force
        shape.torque .+= r × df        # Photon recoil torque
    end
end

function update_thermal_force!(facet::Facet, A_B::Real, ε::Real, n̂, a)
    E = A_B * facet.flux.scat + ε * σ_SB * facet.temps[begin]^4  # Sum of scattering and emission [W/m²]

    @. facet.force = - 2/3 * E * a / c₀ * n̂                      # The first term normal to the face

    for vf in facet.visiblefacets
        @. facet.force += E * a / c₀ * vf.f * vf.d̂               # The second term of self-heating
    end
end



