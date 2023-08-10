
# ****************************************************************
#                     Energy input/output
# ****************************************************************


"""
    total_flux(A_B, A_TH, F_sun, F_scat, F_rad) -> F_total

Total energy absorbed by the face
"""
total_flux(A_B, A_TH, F_sun, F_scat, F_rad) = (1 - A_B) * (F_sun + F_scat) + (1 - A_TH) * F_rad


"""
    energy_io(shape::ShapeModel, params::AbstractThermoParams, nₜ::Integer) -> E_in, E_out, E_cons

Input and output energy per second at a certain time

# Returns
- `E_in`   : Input energy per second at a certain time [W]
- `E_out`  : Output enegey per second at a certain time [W]
- `E_cons` : Output-input energy ratio (`E_out / E_in`)
"""
function energy_io(shape::ShapeModel, params::AbstractThermoParams, nₜ::Integer)
    E_in   = energy_in(shape, params)
    E_out  = energy_out(shape, params, nₜ)
    E_cons = E_out / E_in

    E_in, E_out, E_cons
end


"""
    energy_in(shape::ShapeModel, params::AbstractThermoParams) -> E_in
    energy_in(shape::ShapeModel, A_B)                          -> E_in
    energy_in(A_B, F_sun, F_scat, area)                        -> E_in

Input energy per second on a facet or a whole surface [W]
"""
energy_in(shape::ShapeModel, params::AbstractThermoParams) = energy_in(shape, params.A_B)

function energy_in(shape::ShapeModel, A_B)
    E_in = 0.
    for i in eachindex(shape.faces)
        F_sun = shape.flux[i, 1]
        F_scat = shape.flux[i, 2]
        aᵢ = shape.face_areas[i]

        E_in += energy_in((A_B isa Real ? A_B : A_B[i]), F_sun, F_scat, aᵢ)
    end
    E_in
end

energy_in(A_B, F_sun, F_scat, area) = (1 - A_B) * (F_sun + F_scat) * area


"""
    energy_out(shape::ShapeModel, params::AbstractThermoParams, nₜ::Integer) -> E_out
    energy_out(shape, ε, A_TH, nₜ)                                           -> E_out
    energy_out(ε, T, A_TH, F_rad, area)                                      -> E_out

Output enegey per second from a facet or a whole surface [W]
"""
energy_out(shape::ShapeModel, params::AbstractThermoParams, nₜ::Integer) = energy_out(shape, params.ε, params.A_TH, nₜ)

function energy_out(shape, ε, A_TH, nₜ)
    E_out = 0.
    for i in eachindex(shape.faces)
        Tᵢ = shape.temperature[begin, i, nₜ]
        F_rad = shape.flux[i, 3]
        aᵢ = shape.face_areas[i]

        E_out += energy_out((ε isa Real ? ε : ε[i]), Tᵢ, (A_TH isa Real ? A_TH : A_TH[i]), F_rad, aᵢ)
    end
    E_out
end

energy_out(ε, T, A_TH, F_rad, area) = (ε * σ_SB * T^4 - (1 - A_TH) * F_rad) * area


# ****************************************************************
#                  Energy flux: Sunlight
# ****************************************************************

"""
    update_flux!(shape::ShapeModel, r☉::AbstractVector, params::AbstractThermoParams, nₜ::Integer)

Update energy flux to every facet by solar radiation, scattering, and re-absorption of radiation
"""
function update_flux!(shape::ShapeModel, r☉::AbstractVector, params::AbstractThermoParams, nₜ::Integer)
    update_flux_sun!(shape, r☉)
    update_flux_scat_single!(shape, params)
    update_flux_rad_single!(shape, params, nₜ)
end

"""
    update_flux_sun!(shape, r̂☉, F☉)

Update solar radiation flux on every facet of a shape model.

- `shape` : Shape model
- `r̂☉`    : Normalized vector indicating the direction of the sun in the body-fixed frame
- `F☉`    : Solar radiation flux [W/m²]
"""
function update_flux_sun!(shape::ShapeModel, r̂☉, F☉)
    for i in eachindex(shape.faces)
        if isilluminated(shape, r̂☉, i)
            n̂ = shape.face_normals[i]
            shape.flux[i, 1] = F☉ * (n̂ ⋅ r̂☉)
        else
            shape.flux[i, 1] = 0.
        end
    end
end

"""
    update_flux_sun!(shape, r☉)

Update solar radiation flux on every facet of a shape model.

- `shape` : Shape model
- `r☉`    : Position of the sun in the body-fixed frame, which is not normalized.
"""
function update_flux_sun!(shape::ShapeModel, r☉)
    r̂☉ = SVector{3}(normalize(r☉))
    F☉ = SOLAR_CONST / SPICE.convrt(norm(r☉), "m", "au")^2

    update_flux_sun!(shape, r̂☉, F☉)
end


"""
The secondary is within the critical angle to detect an eclipse event.
"""
function eclipse_is_possible(shapes, sun_from_pri, sec_from_pri)

    R₁ = maximum_radius(shapes[1])
    R₂ = maximum_radius(shapes[2])
    θ_crit = asin((R₁ + R₂) / norm(sec_from_pri))

    r̂☉ = SVector{3}(normalize(sun_from_pri))
    r̂ₛ = SVector{3}(normalize(sec_from_pri))
    θ = acos(r̂☉ ⋅ r̂ₛ)  # Angle of Sun-Primary-Secondary

    θ_crit < θ < π - θ_crit ? false : true
end


"""
    find_eclipse!(shapes, sun_from_pri, sec_from_pri, R₂₁)

- `shapes`
- `sun_from_pri`
- `sec_from_pri` : Position of the secondary relative to primary
- `R₂₁`          : Rotation matrix from secondary to primary
"""
function find_eclipse!(shapes, sun_from_pri, sec_from_pri, R₂₁)

    r̂☉ = SVector{3}(normalize(sun_from_pri))
    # r̂ₛ = SVector{3}(normalize(sec_from_pri))
    rₛ = SVector{3}(sec_from_pri)

    # cosθ = r̂☉ ⋅ r̂ₛ  # Cosine of angle of Sun-Primary-Secondary

    eclipse_is_possible(shapes, r̂☉, rₛ) == false && return

    for i in eachindex(shapes[1].faces)
        shapes[1].flux[i, 1] == 0 && continue  # something wrong?
        A₁, B₁, C₁ = shapes[1].nodes[shapes[1].faces[i]]  # △ABC in primary
        G₁ = shapes[1].face_centers[i]                    # Center of △ABC in primary

        for j in eachindex(shapes[2].faces)
            shapes[2].flux[j, 1] == 0 && continue  # something wrong?
            A₂, B₂, C₂ = shapes[2].nodes[shapes[2].faces[j]]  # △ABC in secondary
            G₂ = shapes[2].face_centers[j]                    # Center of △ABC in secondary
            
            ## Transform coordinates from secondary- to primary-fixed frame
            A₂ = R₂₁ * A₂ + rₛ
            B₂ = R₂₁ * B₂ + rₛ
            C₂ = R₂₁ * C₂ + rₛ
            G₂ = R₂₁ * G₂ + rₛ
            
            raycast(A₂, B₂, C₂, r̂☉, G₁) && (shapes[1].flux[i, 1] = 0)  # something wrong?
            raycast(A₁, B₁, C₁, r̂☉, G₂) && (shapes[2].flux[j, 1] = 0)  # something wrong?
        end
    end
end


# ****************************************************************
#                 Energy flux: Scattering
# ****************************************************************

"""
    update_flux_scat_single!(shape, params::AbstractThermoParams)
    update_flux_scat_single!(shape, A_B)

Update flux of scattered sunlight, only considering single scattering.
"""
update_flux_scat_single!(shape::ShapeModel, params::AbstractThermoParams) = update_flux_scat_single!(shape, params.A_B)

function update_flux_scat_single!(shape, A_B)
    for i in eachindex(shape.faces)
        shape.flux[i, 2] = 0.
        for visiblefacet in shape.visiblefacets[i]
            j = visiblefacet.id
            fᵢⱼ = visiblefacet.f
            A_B = (A_B isa Real ? A_B : A_B[j])
            shape.flux[i, 2] += fᵢⱼ * A_B * shape.flux[j, 1]
        end
    end
end


##= TODO: Implement update_flux_scat_mult! =##

# """
#     update_flux_scat_mult!(shape, params::AbstractThermoParams)
#     update_flux_scat_mult!(shape, A_B)

# Update flux of scattered sunlight, considering multiple scattering.
# """


# ****************************************************************
#                Energy flux: Thermal radiation
# ****************************************************************

"""
    update_flux_rad_single!(shape, params::AbstractThermoParams, nₜ::Integer)
    update_flux_rad_single!(shape, ε, A_TH, nₜ)

Update flux of absorption of thermal radiation from surrounding surface.
Single radiation-absorption is only considered, assuming albedo is close to zero at thermal infrared wavelength.

# Arguments
- `shape`  : Shape model
- `params` : Thermophysical parameters
- `nₜ`     : Index of time step
"""
update_flux_rad_single!(shape::ShapeModel, params::AbstractThermoParams, nₜ::Integer) = update_flux_rad_single!(shape, params.ε, params.A_TH, nₜ)

function update_flux_rad_single!(shape, ε, A_TH, nₜ)
    for i in eachindex(shape.faces)
        shape.flux[i, 3] = 0.
        for visiblefacet in shape.visiblefacets[i]
            j = visiblefacet.id
            fᵢⱼ = visiblefacet.f
            ε = (ε isa Real ? ε : ε[j])
            A_TH = (A_TH isa Real ? A_TH : A_TH[j])
            Tⱼ = shape.temperature[begin, j, nₜ]
            
            shape.flux[i, 3] += ε * σ_SB * (1 - A_TH) * fᵢⱼ * Tⱼ^4
        end
    end
end

