
# ****************************************************************
#                     Energy input/output
# ****************************************************************


"""
    flux_total(A_B, A_TH, F_sun, F_scat, F_rad) -> F_total

Total energy absorbed by the surface [W/m²]

# Arguments
- `A_B`    : Bond albedo [-]
- `A_TH`   : Albedo at thermal infrared wavelength [-]
- `F_sun`  : Flux of direct sunlight [W/m²]
- `F_scat` : Flux of scattered light [W/m²]
- `F_rad`  : Flux of thermal radiation from surrounding surface [W/m²]
"""
flux_total(A_B, A_TH, F_sun, F_scat, F_rad) = (1 - A_B) * F_sun + (1 - A_B) * F_scat + (1 - A_TH) * F_rad


"""
    energy_in(stpm::SingleTPM) -> E_in

Input energy per second on the whole surface [W]
"""
function energy_in(stpm::SingleTPM)
    E_in = 0.
    for nₛ in eachindex(stpm.shape.faces)
        A_B    = (stpm.thermo_params.A_B  isa Real ? stpm.thermo_params.A_B  : stpm.thermo_params.A_B[nₛ])
        A_TH   = (stpm.thermo_params.A_TH isa Real ? stpm.thermo_params.A_TH : stpm.thermo_params.A_TH[nₛ])
        F_sun  = stpm.flux[nₛ, 1]
        F_scat = stpm.flux[nₛ, 2]
        F_rad  = stpm.flux[nₛ, 3]
        a      = stpm.shape.face_areas[nₛ]
        
        E_in += flux_total(A_B, A_TH, F_sun, F_scat, F_rad) * a
    end
    E_in
end


"""
    energy_out(stpm::SingleTPM) -> E_out

Output enegey per second from the whole surface [W]

# Arguments
- `stpm` : Thermophysical model for a single asteroid
"""
function energy_out(stpm::SingleTPM)
    E_out = 0.
    for nₛ in eachindex(stpm.shape.faces)
        ε = (stpm.thermo_params.ε isa Real ? stpm.thermo_params.ε : stpm.thermo_params.ε[nₛ])
        T = stpm.temperature[begin, nₛ]  # Surface temperature
        a = stpm.shape.face_areas[nₛ]

        E_out += ε * σ_SB * T^4 * a
    end
    E_out
end


# ****************************************************************
#                  Energy flux: Sunlight
# ****************************************************************


"""
    update_flux_sun!(stpm::SingleTPM, r̂☉::StaticVector{3}, F☉::Real)

Update solar irradiation flux on every face of a shape model.

- `shape` : Shape model
- `r̂☉`    : Normalized vector indicating the direction of the sun in the body-fixed frame
- `F☉`    : Solar radiation flux [W/m²]
"""
function update_flux_sun!(stpm::SingleTPM, r̂☉::StaticVector{3}, F☉::Real)
    r̂☉ = normalize(r̂☉)

    if stpm.SELF_SHADOWING
        for nₛ in eachindex(stpm.shape.faces)
            if isilluminated(stpm.shape, r̂☉, nₛ)
                n̂ = stpm.shape.face_normals[nₛ]
                stpm.flux[nₛ, 1] = F☉ * (n̂ ⋅ r̂☉)
            else
                stpm.flux[nₛ, 1] = 0
            end
        end
    else
        for nₛ in eachindex(stpm.shape.faces)
            n̂ = stpm.shape.face_normals[nₛ]
            if n̂ ⋅ r̂☉ > 0
                stpm.flux[nₛ, 1] = F☉ * (n̂ ⋅ r̂☉)
            else
                stpm.flux[nₛ, 1] = 0
            end
        end
    end
end


"""
    update_flux_sun!(stpm::SingleTPM, r☉::StaticVector{3})

Update solar irradiation flux on every face of a shape model.

# Arguments
- `stpm` : Thermophysical model for a single asteroid
- `r☉`   : Position of the sun in the body-fixed frame (NOT normalized)
"""
function update_flux_sun!(stpm::SingleTPM, r☉::StaticVector{3})
    r̂☉ = SVector{3}(normalize(r☉))
    F☉ = SOLAR_CONST / SPICE.convrt(norm(r☉), "m", "au")^2

    update_flux_sun!(stpm, r̂☉, F☉)
end


"""
    update_flux_sun!(btpm::BinaryTPM, r☉₁::StaticVector{3}, r☉₂::StaticVector{3})

# Arguments
- `btpm` : Thermophysical model for a binary asteroid
- `r☉₁`  : Sun's position in the body-fixed frame of the primary, which is not normalized.
- `r☉₂`  : Sun's position in the body-fixed frame of the secondary, which is not normalized.
"""
function update_flux_sun!(btpm::BinaryTPM, r☉₁::StaticVector{3}, r☉₂::StaticVector{3})
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
    stpm.SELF_HEATING == false && return

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
    update_flux_rad_single!(stpm::SingleTPM)

Update flux of absorption of thermal radiation from surrounding surface.
Single radiation-absorption is only considered, assuming albedo is close to zero at thermal infrared wavelength.

# Arguments
- `stpm` : Thermophysical model for a single asteroid
"""
function update_flux_rad_single!(stpm::SingleTPM)
    stpm.SELF_HEATING == false && return

    for nₛ in eachindex(stpm.shape.faces)
        stpm.flux[nₛ, 3] = 0.
        for visiblefacet in stpm.shape.visiblefacets[nₛ]
            j    = visiblefacet.id
            fᵢⱼ  = visiblefacet.f
            ε    = (stpm.thermo_params.ε    isa Real ? stpm.thermo_params.ε    : stpm.thermo_params.ε[j])
            A_TH = (stpm.thermo_params.A_TH isa Real ? stpm.thermo_params.A_TH : stpm.thermo_params.A_TH[j])
            Tⱼ   = stpm.temperature[begin, j]
            
            stpm.flux[nₛ, 3] += ε * σ_SB * (1 - A_TH) * fᵢⱼ * Tⱼ^4
        end
    end
end


"""
    update_flux_rad_single!(btpm::BinaryTPM)

Update flux of absorption of thermal radiation from surrounding surface.
Single radiation-absorption is only considered, assuming albedo is close to zero at thermal infrared wavelength.

# Arguments
- `btpm` : Thermophysical model for a binary asteroid
"""
function update_flux_rad_single!(btpm::BinaryTPM)
    update_flux_rad_single!(btpm.pri)
    update_flux_rad_single!(btpm.sec)
end


# ****************************************************************
#                Mutual shadowing of binary asteroid
# ****************************************************************


"""
    mutual_shadowing!(btpm::BinaryTPM, r☉, rₛ, R₂₁)

Detect eclipse events between the primary and secondary, and update the solar fluxes of the faces.

# Arguments
- `btpm` : Thermophysical model for a binary asteroid
- `r☉`   : Position of the sun relative to the primary       (NOT normalized)
- `rₛ`   : Position of the secondary relative to the primary (NOT normalized)
- `R₂₁`  : Rotation matrix from secondary to primary
"""
function mutual_shadowing!(btpm::BinaryTPM, r☉, rₛ, R₂₁)
    btpm.MUTUAL_SHADOWING == false && return

    shape1 = btpm.pri.shape
    shape2 = btpm.sec.shape
    r̂☉ = normalize(r☉)

    cosθ = normalize(r☉) ⋅ normalize(rₛ)  # Cosine of the angle of Sun-Primary-Secondary
    cosθ = min(1.0, max(-1.0, cosθ))      # Clip the value to [-1.0, 1.0]
    θ = acos(cosθ)

    R₁ = maximum_radius(shape1)
    R₂ = maximum_radius(shape2)
    R₁ < R₂ && error("Error: The primary radius is smaller than the secondary.")
    
    sinθ₊ = (R₁ + R₂) / norm(rₛ)        # Sine of the critical angle at which partial ecripse can occur
    sinθ₋ = (R₁ - R₂) / norm(rₛ)        # Sine of the critical angle at which total ecripse can occur
    sinθ₊ = min(1.0, max(-1.0, sinθ₊))  # Clip the value to [-1.0, 1.0]
    sinθ₋ = min(1.0, max(-1.0, sinθ₋))  # Clip the value to [-1.0, 1.0]

    θ₊ = asin(sinθ₊)  # Critical angle at which partial ecripse can occur
    θ₋ = asin(sinθ₋)  # Critical angle at which total ecripse can occur

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

            d₁ₛ = rₛ - G₁                                        # Vector from △A₁B₁C₁ to secondary center
            θ₁ = acos(min(1.0, max(-1.0, r̂☉ ⋅ normalize(d₁ₛ))))  # Angle of Sun-△A₁B₁C₁-Secondary
            θ_R₂ = asin(min(1.0, max(-1.0, R₂ / norm(d₁ₛ))))     # Critical angle related to the maximum radius of the secondary
            θ_r₂ = asin(min(1.0, max(-1.0, r₂ / norm(d₁ₛ))))     # Critical angle related to the minimum radius of the secondary

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

            d₂ₚ = - G₂                                           # Vector from △A₂B₂C₂ to primary center (origin)
            θ₂ = acos(min(1.0, max(-1.0, r̂☉ ⋅ normalize(d₂ₚ))))  # Angle of Sun-△A₂B₂C₂-Primary
            θ_R₁ = asin(min(1.0, max(-1.0, R₁ / norm(d₂ₚ))))     # Critical angle related to the maximum radius of the primary
            θ_r₁ = asin(min(1.0, max(-1.0, r₁ / norm(d₂ₚ))))     # Critical angle related to the minimum radius of the primary

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


# ****************************************************************
#                Mutual heating of binary asteroid
# ****************************************************************


"""
    mutual_heating!(btpm::BinaryTPM, rₛ, R₂₁)

Calculate the mutual heating between the primary and secondary asteroids.

# Arguments
- `btpm` : Thermophysical model for a binary asteroid
- `rₛ`   : Position of the secondary relative to the primary (NOT normalized)
- `R₂₁`  : Rotation matrix from secondary to primary

# TO DO
- Need to consider local horizon?
"""
function mutual_heating!(btpm::BinaryTPM, rₛ, R₂₁)
    btpm.MUTUAL_HEATING == false && return

    shape1 = btpm.pri.shape
    shape2 = btpm.sec.shape
    thermo_params1 = btpm.pri.thermo_params
    thermo_params2 = btpm.sec.thermo_params

    for i in eachindex(shape1.faces)  # △A₁B₁C₁ in primary
        c₁ = shape1.face_centers[i]   # Center of △A₁B₁C₁
        n̂₁ = shape1.face_normals[i]   # Normal vector of △A₁B₁C₁
        a₁ = shape1.face_areas[i]     # Area of △A₁B₁C₁

        for j in eachindex(shape2.faces)  # △A₂B₂C₂ in secondary
            c₂ = shape2.face_centers[j]   # Center of △A₂B₂C₂
            n̂₂ = shape2.face_normals[j]   # Normal vector of △A₂B₂C₂
            a₂ = shape2.face_areas[j]     # Area of △A₂B₂C₂
        
            ## Transformation from secondary to primary frame
            c₂ = R₂₁ * c₂ + rₛ
            n̂₂ = R₂₁ * n̂₂

            f₁₂, d₁₂, d̂₁₂ = view_factor(c₁, c₂, n̂₁, n̂₂, a₂)  # View factor from △A₁B₁C₁ to △A₂B₂C₂
            f₂₁, d₂₁, d̂₂₁ = view_factor(c₂, c₁, n̂₂, n̂₁, a₁)  # View factor from △A₂B₂C₂ to △A₁B₁C₁

            ## if △A₁B₁C₁ and △A₂B₂C₂ are facing each other
            if d̂₁₂ ⋅ n̂₁ > 0 && d̂₁₂ ⋅ n̂₂ < 0
                T₁ = btpm.pri.temperature[begin, i]
                T₂ = btpm.sec.temperature[begin, j]

                ε₁    = (thermo_params1.ε    isa Real ? thermo_params1.ε    : thermo_params1.ε[i])
                ε₂    = (thermo_params2.ε    isa Real ? thermo_params2.ε    : thermo_params2.ε[j])
                A_B₁  = (thermo_params1.A_B  isa Real ? thermo_params1.A_B  : thermo_params1.A_B[i])
                A_B₂  = (thermo_params2.A_B  isa Real ? thermo_params2.A_B  : thermo_params2.A_B[j])
                A_TH₁ = (thermo_params1.A_TH isa Real ? thermo_params1.A_TH : thermo_params1.A_TH[i])
                A_TH₂ = (thermo_params2.A_TH isa Real ? thermo_params2.A_TH : thermo_params2.A_TH[j])

                ## Mutual heating by scattered light
                btpm.pri.flux[i, 2] += f₁₂ * A_B₂ * btpm.sec.flux[j, 1]
                btpm.sec.flux[j, 2] += f₂₁ * A_B₁ * btpm.pri.flux[i, 1]

                ## Mutual heating by thermal radiation
                btpm.pri.flux[i, 3] += ε₂ * σ_SB * (1 - A_TH₂) * f₁₂ * T₂^4
                btpm.sec.flux[j, 3] += ε₁ * σ_SB * (1 - A_TH₁) * f₂₁ * T₁^4
            end
        end
    end
end
