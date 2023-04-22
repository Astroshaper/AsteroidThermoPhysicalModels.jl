

# ****************************************************************
#                  Photon recoil force / torque
# ****************************************************************

"""
    update_force!(shape::ShapeModel, params)
    update_force!(shape::ShapeModel, A_B::Real,           ϵ::Real)
    update_force!(shape::ShapeModel, A_B::AbstractVector, ϵ::Real)
    update_force!(shape::ShapeModel, A_B::Real,           ϵ::AbstractVector)
    update_force!(shape::ShapeModel, A_B::AbstractVector, ϵ::AbstractVector)
    update_force!(facet::Facet,      A_B::Real,           ϵ::Real)

Update photon recoil force on every facet (dfᵢ)
"""
update_force!(shape::ShapeModel, params) = update_force!(shape, params.A_B, params.ϵ)

function update_force!(shape::ShapeModel, A_B::Real, ϵ::Real)
    for facet in shape.facets
        update_force!(facet, A_B, ϵ)
    end
end

function update_force!(shape::ShapeModel, A_B::AbstractVector, ϵ::Real)
    for (i, facet) in enumerate(shape.facets)
        update_force!(facet, A_B[i], ϵ)
    end
end

function update_force!(shape::ShapeModel, A_B::Real, ϵ::AbstractVector)
    for (i, facet) in enumerate(shape.facets)
        update_force!(facet, A_B, ϵ[i])
    end
end

function update_force!(shape::ShapeModel, A_B::AbstractVector, ϵ::AbstractVector)
    for (i, facet) in enumerate(shape.facets)
        update_force!(facet, A_B[i], ϵ[i])
    end
end

function update_force!(facet::Facet, A_B::Real, ϵ::Real)
    E = A_B * facet.flux.scat + ϵ * σ_SB * facet.temps[begin]^4

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

"""
    update_force_Rubincam!(shape, params_thermo)

Update photon recoil force on every facet (df) based on Rubincam (2000) approximation
"""
function update_force_Rubincam!(shape, params_thermo)
    shape.torque .= 0.
    for facet in shape.facets
        facet.force .= - 2 * facet.flux.sun * facet.area / (3*c₀) .* facet.normal
        shape.torque .+= facet.center × SVector{3}(facet.force)
    end
end