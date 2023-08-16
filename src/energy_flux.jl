
# ****************************************************************
#                     Energy input/output
# ****************************************************************


"""
    flux_total(A_B, A_TH, F_sun, F_scat, F_rad) -> F_total

Total energy absorbed by the surface

# Arguments
- `A_B`    : Bond albedo
- `A_TH`   : Albedo at thermal infrared wavelength
- `F_sun`  : Flux of direct sunlight
- `F_scat` : Flux of scattered light
- `F_rad`  : Flux of thermal radiation from surrounding surface
"""
flux_total(A_B, A_TH, F_sun, F_scat, F_rad) = (1 - A_B) * (F_sun + F_scat) + (1 - A_TH) * F_rad


"""
    energy_io(stpm::SingleTPM, nₜ::Integer) -> E_in, E_out, E_cons

Input and output energy per second at a given time step

# Returns
- `E_in`   : Input energy per second at a given time step [W]
- `E_out`  : Output enegey per second at a given time step [W]
- `E_cons` : Output-input energy ratio (`E_out / E_in`)
"""
function energy_io(stpm::SingleTPM, nₜ::Integer)
    E_in   = energy_in(stpm)
    E_out  = energy_out(stpm, nₜ)
    E_cons = E_out / E_in

    E_in, E_out, E_cons
end


"""
    energy_in(stpm::SingleTPM) -> E_in

Input energy per second on a facet or a whole surface [W]
"""
energy_in(shape::ShapeModel, params::AbstractThermoParams) = energy_in(shape, params.A_B)

function energy_in(stpm::SingleTPM)
    E_in = 0.
    for i in eachindex(stpm.shape.faces)
        A_B    = (stpm.thermo_params.A_B isa Real ? stpm.thermo_params.A_B : stpm.thermo_params.A_B[i])
        F_sun  = stpm.flux[i, 1]
        F_scat = stpm.flux[i, 2]
        aᵢ     = stpm.shape.face_areas[i]
        
        E_in += energy_in(A_B, F_sun, F_scat, aᵢ)
    end
    E_in
end


"""
    energy_in(A_B, F_sun, F_scat, area) -> E_in

Input energy per second on a face [W]

# Arguments
- `A_B`    : Bond albedo
- `F_sun`  : Flux of direct sunlight
- `F_scat` : Flux of scattered light
- `area`   : Area of the face
"""
energy_in(A_B, F_sun, F_scat, area) = (1 - A_B) * (F_sun + F_scat) * area


"""
    energy_out(stpm::SingleTPM, nₜ::Integer) -> E_out

Output enegey per second from a entire surface [W]

# Arguments
- `stpm` : Thermophysical model for a single asteroid
- `nₜ`   : Index of time step
"""
function energy_out(stpm::SingleTPM, nₜ::Integer)
    E_out = 0.
    for i in eachindex(stpm.shape.faces)
        Tᵢ    = stpm.temperature[begin, i, nₜ]
        F_rad = stpm.flux[i, 3]
        aᵢ    = stpm.shape.face_areas[i]
        ε     = (stpm.thermo_params.ε    isa Real ? stpm.thermo_params.ε    : stpm.thermo_params.ε[i])
        A_TH  = (stpm.thermo_params.A_TH isa Real ? stpm.thermo_params.A_TH : stpm.thermo_params.A_TH[i])

        E_out += energy_out(ε, Tᵢ, A_TH, F_rad, aᵢ)
    end
    E_out
end


"""
    energy_out(ε, T, A_TH, F_rad, area) -> E_out

Output enegey per second from a face [W]

# Arguments
- `ε`     : Emissivity
- `T`     : Temperature of the face
- `A_TH`  : Albedo at thermal infrared wavelength
- `F_rad` : Flux of thermal radiation from surrounding surface
- `area`  : Area of the face
"""
energy_out(ε, T, A_TH, F_rad, area) = (ε * σ_SB * T^4 - (1 - A_TH) * F_rad) * area


# ****************************************************************
#                  Energy flux: Sunlight
# ****************************************************************


"""
    update_flux_sun!(stpm::SingleTPM, r̂☉::AbstractVector, F☉::Real)

Update solar irradiation flux on every face of a shape model.

- `shape` : Shape model
- `r̂☉`    : Normalized vector indicating the direction of the sun in the body-fixed frame
- `F☉`    : Solar radiation flux [W/m²]
"""
function update_flux_sun!(stpm::SingleTPM, r̂☉::AbstractVector, F☉::Real)
    r̂☉ = SVector{3}(normalize(r̂☉))
    for nₛ in eachindex(stpm.shape.faces)
        if isilluminated(stpm.shape, r̂☉, nₛ)
            n̂ = stpm.shape.face_normals[nₛ]
            stpm.flux[nₛ, 1] = F☉ * (n̂ ⋅ r̂☉)
        else
            stpm.flux[nₛ, 1] = 0.
        end
    end
end


"""
    update_flux_sun!(stpm::SingleTPM, r☉::AbstractVector)

Update solar irradiation flux on every face of a shape model.

# Arguments
- `stpm` : Thermophysical model for a single asteroid
- `r☉`   : Position of the sun in the body-fixed frame, which is not normalized.
"""
function update_flux_sun!(stpm::SingleTPM, r☉::AbstractVector)
    r̂☉ = SVector{3}(normalize(r☉))
    F☉ = SOLAR_CONST / SPICE.convrt(norm(r☉), "m", "au")^2

    update_flux_sun!(stpm, r̂☉, F☉)
end


"""
    update_flux_sun!(btpm::BinaryTPM, r☉₁::AbstractVector, r☉₂::AbstractVector)

# Arguments
- `btpm` : Thermophysical model for a binary asteroid
- `r☉₁`  : Sun's position in the body-fixed frame of the primary, which is not normalized.
- `r☉₂`  : Sun's position in the body-fixed frame of the secondary, which is not normalized.
"""
function update_flux_sun!(btpm::BinaryTPM, r☉₁::AbstractVector, r☉₂::AbstractVector)
    update_flux_sun!(btpm.pri, r☉₁)
    update_flux_sun!(btpm.sec, r☉₂)
end


# ****************************************************************
#                 Energy flux: Scattering
# ****************************************************************

"""
    update_flux_scat_single!(stpm::SingleTPM)

Update flux of scattered sunlight, only considering single scattering.

# Arguments
- `stpm` : Thermophysical model for a single asteroid
"""
function update_flux_scat_single!(stpm::SingleTPM)
    for nₛ in eachindex(stpm.shape.faces)
        stpm.flux[nₛ, 2] = 0.
        for visiblefacet in stpm.shape.visiblefacets[nₛ]
            j   = visiblefacet.id
            fᵢⱼ = visiblefacet.f
            A_B = (stpm.thermo_params.A_B isa Real ? stpm.thermo_params.A_B : stpm.thermo_params.A_B[j])

            stpm.flux[nₛ, 2] += fᵢⱼ * A_B * stpm.flux[j, 1]
        end
    end
end


"""
    update_flux_scat_single!(btpm::BinaryTPM)

Update flux of scattered sunlight, only considering single scattering.

# Arguments
- `btpm` : Thermophysical model for a binary asteroid
"""
function update_flux_scat_single!(btpm::BinaryTPM)
    update_flux_scat_single!(btpm.pri)
    update_flux_scat_single!(btpm.sec)
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
    update_flux_rad_single!(stpm::SingleTPM, nₜ::Integer)

Update flux of absorption of thermal radiation from surrounding surface.
Single radiation-absorption is only considered, assuming albedo is close to zero at thermal infrared wavelength.

# Arguments
- `stpm` : Thermophysical model for a single asteroid
- `nₜ`   : Index of time step
"""
function update_flux_rad_single!(stpm::SingleTPM, nₜ::Integer)
    for nₛ in eachindex(stpm.shape.faces)
        stpm.flux[nₛ, 3] = 0.
        for visiblefacet in stpm.shape.visiblefacets[nₛ]
            j    = visiblefacet.id
            fᵢⱼ  = visiblefacet.f
            ε    = (stpm.thermo_params.ε    isa Real ? stpm.thermo_params.ε    : stpm.thermo_params.ε[j])
            A_TH = (stpm.thermo_params.A_TH isa Real ? stpm.thermo_params.A_TH : stpm.thermo_params.A_TH[j])
            Tⱼ   = stpm.temperature[begin, j, nₜ]
            
            stpm.flux[nₛ, 3] += ε * σ_SB * (1 - A_TH) * fᵢⱼ * Tⱼ^4
        end
    end
end


"""
    update_flux_rad_single!(btpm::BinaryTPM, nₜ::Integer)

Update flux of absorption of thermal radiation from surrounding surface.
Single radiation-absorption is only considered, assuming albedo is close to zero at thermal infrared wavelength.

# Arguments
- `btpm` : Thermophysical model for a binary asteroid
- `nₜ`   : Index of time step
"""
function update_flux_rad_single!(btpm::BinaryTPM, nₜ::Integer)
    update_flux_rad_single!(btpm.pri, nₜ)
    update_flux_rad_single!(btpm.sec, nₜ)
end


# ****************************************************************
#                  Eclipse of binary asteroid
# ****************************************************************

"""
The secondary is within the critical angle to detect an eclipse event.
"""
function binary_is_aligned(btpm::BinaryTPM, sun_from_pri, sec_from_pri)

    R₁ = maximum_radius(btpm.pri.shape)
    R₂ = maximum_radius(btpm.sec.shape)
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
function find_eclipse!(btpm::BinaryTPM, sun_from_pri, sec_from_pri, R₂₁)

    r̂☉ = SVector{3}(normalize(sun_from_pri))
    # r̂ₛ = SVector{3}(normalize(sec_from_pri))
    rₛ = SVector{3}(sec_from_pri)

    binary_is_aligned(btpm, r̂☉, rₛ) == false && return

    for i in eachindex(btpm.pri.shape.faces)
        btpm.pri.flux[i, 1] == 0 && continue                        # something wrong?
        A₁, B₁, C₁ = btpm.pri.shape.nodes[btpm.pri.shape.faces[i]]  # △ABC in primary
        G₁ = btpm.pri.shape.face_centers[i]                         # Center of △ABC in primary

        for j in eachindex(btpm.sec.shape.faces)
            btpm.sec.flux[j, 1] == 0 && continue                        # something wrong?
            A₂, B₂, C₂ = btpm.sec.shape.nodes[btpm.sec.shape.faces[j]]  # △ABC in secondary
            G₂ = btpm.sec.shape.face_centers[j]                         # Center of △ABC in secondary
            
            ## Transform coordinates from secondary- to primary-fixed frame
            A₂ = R₂₁ * A₂ + rₛ
            B₂ = R₂₁ * B₂ + rₛ
            C₂ = R₂₁ * C₂ + rₛ
            G₂ = R₂₁ * G₂ + rₛ
            
            raycast(A₂, B₂, C₂, r̂☉, G₁) && (btpm.pri.flux[i, 1] = 0)  # something wrong?
            raycast(A₁, B₁, C₁, r̂☉, G₂) && (btpm.sec.flux[j, 1] = 0)  # something wrong?
        end
    end
end
