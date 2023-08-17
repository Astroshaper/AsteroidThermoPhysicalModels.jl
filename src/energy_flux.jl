
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
    find_eclipse!(btpm::BinaryTPM, r☉, rₛ, R₂₁)

Find eclipse events between the primary and secondary, and update the solar fluxes of the faces.

# Arguments
- `btpm` : Thermophysical model for a binary asteroid
- `r☉`   : Position of the sun relative to the primary       (NOT normalized)
- `rₛ`   : Position of the secondary relative to the primary (NOT normalized)
- `R₂₁`  : Rotation matrix from secondary to primary

# TO DO:
- r₁ = minimum_radius(shape1)
- r₂ = minimum_radius(shape2)
Use these radii and conduct early out of shadowed faces before calling `raycast`.
"""
function find_eclipse!(btpm::BinaryTPM, r☉, rₛ, R₂₁)

    shape1 = btpm.pri.shape
    shape2 = btpm.sec.shape
    r̂☉ = normalize(r☉)

    θ = acos((r☉ ⋅ rₛ) / (norm(r☉) * norm(rₛ)))  # Angle of Sun-Primary-Secondary

    R₁ = maximum_radius(shape1)
    R₂ = maximum_radius(shape2)
    R₁ < R₂ && error("Error: The primary radius is smaller than the secondary.")
    
    θ₊ = asin((R₁ + R₂) / norm(rₛ))  # Critical angle at which partial ecripse can occur
    θ₋ = asin((R₁ - R₂) / norm(rₛ))  # Critical angle at which total ecripse can occur

    #### Partital eclipse of the primary ####
    if 0 ≤ θ < θ₊
        r₂ = minimum_radius(shape2)

        for i in eachindex(shape1.faces)
            G₁ = shape1.face_centers[i]  # Center of △A₁B₁C₁ in primary
            n̂₁ = shape1.face_normals[i]  # Normal vector of △A₁B₁C₁ in primary

            ## if △A₁B₁C₁ is NOT facing the sun
            if r̂☉ ⋅ n̂₁ < 0
                btpm.pri.flux[i, 1] = 0
                continue
            end

            d₁ₛ = rₛ - G₁                   # Vector from △A₁B₁C₁ to secondary center
            θ₁ = acos(r̂☉ ⋅ normalize(d₁ₛ))  # Angle of Sun-△A₁B₁C₁-Secondary
            θ_R₂ = asin(R₂ / norm(d₁ₛ))     # Critical angle related to the maximum radius of the secondary
            θ_r₂ = asin(r₂ / norm(d₁ₛ))     # Critical angle related to the minimum radius of the secondary

            ## In the secondary shadow
            if θ₁ < θ_r₂
                btpm.pri.flux[i, 1] = 0
                continue
            ## Out of the secondary shadow
            elseif θ₁ > θ_R₂
                continue
            else
                for j in eachindex(shape2.faces)
                    A₂, B₂, C₂ = shape2.nodes[shape2.faces[j]]  # △A₂B₂C₂ in secondary
                    G₂ = shape2.face_centers[j]                 # Center of △A₂B₂C₂
                    n̂₂ = shape2.face_normals[j]                 # Normal vector of △A₂B₂C₂
            
                    ## Transformation from secondary to primary frame
                    A₂ = R₂₁ * A₂ + rₛ
                    B₂ = R₂₁ * B₂ + rₛ
                    C₂ = R₂₁ * C₂ + rₛ
                    G₂ = R₂₁ * G₂ + rₛ
                    n̂₂ = R₂₁ * n̂₂

                    d₁₂ = G₂ - G₁  # Vector from primary face i to secondary face j
                
                    ## if △A₁B₁C₁ and △A₂B₂C₂ are facing each other
                    if d₁₂ ⋅ n̂₁ > 0 && d₁₂ ⋅ n̂₂ < 0
                        if raycast(A₂, B₂, C₂, r̂☉, G₁)
                            btpm.pri.flux[i, 1] = 0
                            break
                        end
                    end
                end
            end
        end
    
    #### No eclipse ####
    # elseif θ₊ ≤ θ < π - θ₊
    # Do nothing
    
    #### Partial eclipse of the secondary ####
    elseif π - θ₊ ≤ θ < π - θ₋
        r₁ = minimum_radius(shape1)

        for j in eachindex(shape2.faces)
            G₂ = shape2.face_centers[j]  # Center of △A₂B₂C₂ in secondary
            n̂₂ = shape2.face_normals[j]  # Normal vector of △A₂B₂C₂ in secondary
        
            ## Transformation from secondary to primary frame
            G₂ = R₂₁ * G₂ + rₛ
            n̂₂ = R₂₁ * n̂₂

            ## if △A₂B₂C₂ is NOT facing the sun
            if r̂☉ ⋅ n̂₂ < 0
                btpm.sec.flux[j, 1] = 0
                continue
            end

            d₂ₚ = - G₂                      # Vector from △A₂B₂C₂ to primary center (origin)
            θ₂ = acos(r̂☉ ⋅ normalize(d₂ₚ))  # Angle of Sun-△A₂B₂C₂-Primary
            θ_R₁ = asin(R₁ / norm(d₂ₚ))     # Critical angle related to the maximum radius of the primary
            θ_r₁ = asin(r₁ / norm(d₂ₚ))     # Critical angle related to the minimum radius of the primary

            ## In the primary shadow
            if θ₂ < θ_r₁
                btpm.sec.flux[j, 1] = 0
                continue
            ## Out of the primary shadow
            elseif θ₂ > θ_R₁
                continue
            else
                for i in eachindex(shape1.faces)
                    A₁, B₁, C₁ = shape1.nodes[shape1.faces[i]]  # △A₁B₁C₁ in primary
                    G₁ = shape1.face_centers[i]                 # Center of △A₁B₁C₁
                    n̂₁ = shape1.face_normals[i]                 # Normal vector of △A₁B₁C₁

                    d₁₂ = G₂ - G₁  # Vector from primary face i to secondary face j
                
                    ## if △A₁B₁C₁ and △A₂B₂C₂ are facing each other
                    if d₁₂ ⋅ n̂₁ > 0 && d₁₂ ⋅ n̂₂ < 0
                        if raycast(A₁, B₁, C₁, r̂☉, G₂)
                            btpm.sec.flux[j, 1] = 0
                            break
                        end
                    end
                end
            end
        end
    
    #### Total eclipse of the secondary ####
    elseif π - θ₋ ≤ θ < π
        btpm.sec.flux[:, 1] .= 0
    end
end

