

# ****************************************************************
#                  Photon recoil force / torque
# ****************************************************************

"""
    update_force!(shape::ShapeModel, params)
    update_force!(shape::ShapeModel, A_B::Real,           ε::Real)
    update_force!(shape::ShapeModel, A_B::AbstractVector, ε::Real)
    update_force!(shape::ShapeModel, A_B::Real,           ε::AbstractVector)
    update_force!(shape::ShapeModel, A_B::AbstractVector, ε::AbstractVector)
    update_force!(facet::Facet,      A_B::Real,           ε::Real)

Update photon recoil force on every facet (dfᵢ)
"""
update_force!(shape::ShapeModel, params) = update_force!(shape, params.A_B, params.ε)

function update_force!(shape::ShapeModel, A_B::Real, ε::Real)
    for facet in shape.facets
        update_force!(facet, A_B, ε)
    end
end

function update_force!(shape::ShapeModel, A_B::AbstractVector, ε::Real)
    for (i, facet) in enumerate(shape.facets)
        update_force!(facet, A_B[i], ε)
    end
end

function update_force!(shape::ShapeModel, A_B::Real, ε::AbstractVector)
    for (i, facet) in enumerate(shape.facets)
        update_force!(facet, A_B, ε[i])
    end
end

function update_force!(shape::ShapeModel, A_B::AbstractVector, ε::AbstractVector)
    for (i, facet) in enumerate(shape.facets)
        update_force!(facet, A_B[i], ε[i])
    end
end

function update_force!(facet::Facet, A_B::Real, ε::Real)
    E = A_B * facet.flux.scat + ε * σ_SB * facet.temps[begin]^4

    @. facet.force = facet.normal
    for vf in facet.visiblefacets
        @. facet.force -= 3/2 * vf.f * vf.d̂
    end
    @. facet.force *= - 2/3 * E * facet.area / c₀
end

"""
    sum_force_torque!(shape::Shape)

Integrate the force and torque over the global surface
"""
function sum_force_torque!(shape::ShapeModel)
    shape.force  .= 0
    shape.torque .= 0
    for facet in shape.facets
        r  = SVector{3}(facet.center)
        r̂  = normalize(r)
        df = SVector{3}(facet.force)
        
        shape.force  .+= (r̂ ⋅ df) * r̂  # Photon recoil force
        shape.torque .+= r × df        # Photon recoil torque
    end
end

