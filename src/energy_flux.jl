#=
energy_flux.jl

Energy flux calculations for asteroid thermophysical modeling.
This file contains functions for computing various energy fluxes including:
- Direct solar radiation
- Scattered sunlight from other faces
- Thermal radiation from surrounding surfaces
- Total energy input/output balance
- Mutual shadowing and heating effects for binary asteroids
=#

# ╔═══════════════════════════════════════════════════════════════════╗
# ║                     Energy input/output                           ║
# ╚═══════════════════════════════════════════════════════════════════╝

"""
    absorbed_energy_flux(R_vis, R_ir, F_sun, F_scat, F_rad) -> F_abs

Calculate the total energy flux absorbed by a surface element, accounting for
wavelength-dependent reflectance properties.

# Arguments
- `R_vis::Real` : Reflectance for visible light [-], valid between 0 and 1.
- `R_ir::Real` : Reflectance for thermal infrared [-], valid between 0 and 1.
- `F_sun::Real` : Direct solar radiation flux [W/m²]
- `F_scat::Real` : Scattered sunlight flux from other surfaces [W/m²]
- `F_rad::Real` : Thermal radiation flux from surrounding surfaces [W/m²]

# Returns
- `F_abs::Real` : Total absorbed energy flux [W/m²]

# Mathematical Formula
```
F_abs = (1 - R_vis) × F_sun + (1 - R_vis) × F_scat + (1 - R_ir) × F_rad
```

# Physical Interpretation
The function accounts for different reflectance properties at different wavelengths:
- Solar radiation (F_sun) and scattered light (F_scat) are in the visible spectrum
- Thermal radiation (F_rad) is in the infrared spectrum
- The absorbed fraction is (1 - reflectance) for each component

# Example
```julia
R_vis = 0.1   # 10% reflectance in visible
R_ir = 0.05   # 5% reflectance in IR
F_sun = 1000.0   # Direct solar flux
F_scat = 50.0    # Scattered light
F_rad = 100.0    # Thermal radiation
F_abs = absorbed_energy_flux(R_vis, R_ir, F_sun, F_scat, F_rad)
# Returns: 0.9 × 1000 + 0.9 × 50 + 0.95 × 100 = 1040.0 W/m²
```
"""
absorbed_energy_flux(R_vis, R_ir, F_sun, F_scat, F_rad) = (1 - R_vis) * F_sun + (1 - R_vis) * F_scat + (1 - R_ir) * F_rad


"""
    energy_in(stpm::SingleAsteroidTPM) -> E_in

Calculate the total energy input rate (power) absorbed by the entire asteroid surface.

# Arguments
- `stpm::SingleAsteroidTPM` : Thermophysical model for a single asteroid

# Returns
- `E_in::Float64` : Total absorbed power [W]

# Calculation
Integrates the absorbed energy flux over all surface facets:
```
E_in = Σᵢ F_abs,ᵢ × Aᵢ
```
where:
- F_abs,ᵢ is the absorbed energy flux on facet i (calculated by `absorbed_energy_flux`)
- Aᵢ is the area of facet i

# Components
The absorbed energy includes:
1. Direct solar radiation
2. Scattered light from other facets (if self-heating is enabled)
3. Thermal radiation from other facets (if self-heating is enabled)

# Usage
This function is typically used to check energy conservation in the model by
comparing with `energy_out`.

# See Also
- `energy_out` for the total emitted power
- `absorbed_energy_flux` for the flux calculation
"""
function energy_in(stpm::SingleAsteroidTPM)
    E_in = 0.
    for i in eachindex(stpm.shape.faces)
        R_vis  = stpm.thermo_params.reflectance_vis[i]
        R_ir   = stpm.thermo_params.reflectance_ir[i]
        F_sun  = stpm.flux_sun[i]
        F_scat = stpm.flux_scat[i]
        F_rad  = stpm.flux_rad[i]
        a      = stpm.shape.face_areas[i]
        
        E_in += absorbed_energy_flux(R_vis, R_ir, F_sun, F_scat, F_rad) * a
    end
    E_in
end


"""
    energy_out(stpm::SingleAsteroidTPM) -> E_out

Calculate the total energy output rate (power) emitted by the entire asteroid surface
through thermal radiation.

# Arguments
- `stpm::SingleAsteroidTPM` : Thermophysical model for a single asteroid

# Returns
- `E_out::Float64` : Total emitted power [W]

# Calculation
Integrates the thermal emission over all surface facets using the Stefan-Boltzmann law:
```
E_out = Σᵢ εᵢ × σ × Tᵢ⁴ × Aᵢ
```
where:
- εᵢ is the emissivity of facet i
- σ is the Stefan-Boltzmann constant (5.67×10⁻⁸ W/m²/K⁴)
- Tᵢ is the surface temperature of facet i [K]
- Aᵢ is the area of facet i [m²]

# Energy Conservation
In thermal equilibrium, E_out should approximately equal E_in (from `energy_in`).
The ratio E_out/E_in is often used as a convergence criterion in thermophysical models.

# See Also
- `energy_in` for the total absorbed power
- `update_thermal_force!` for thermal recoil effects from this emission
"""
function energy_out(stpm::SingleAsteroidTPM)
    E_out = 0.
    for i in eachindex(stpm.shape.faces)
        ε = stpm.thermo_params.emissivity[i]
        T = stpm.temperature[begin, i]  # Surface temperature
        a = stpm.shape.face_areas[i]

        E_out += ε * σ_SB * T^4 * a
    end
    E_out
end


# ╔═══════════════════════════════════════════════════════════════════╗
# ║                     Unified flux update API                       ║
# ╚═══════════════════════════════════════════════════════════════════╝

"""
    update_flux_all!(stpm::SingleAsteroidTPM, r☉::StaticVector{3})

Update all energy fluxes (solar, scattered, thermal radiation) to the surface for a single asteroid.

# Arguments
- `stpm::SingleAsteroidTPM` : Thermophysical model for a single asteroid
- `r☉::StaticVector{3}`     : Sun's position in the asteroid-fixed frame (NOT normalized) [m]

# Algorithm
1. Updates direct solar flux on all faces considering self-shadowing
2. Updates scattered sunlight flux from other faces (self-heating)
3. Updates thermal radiation flux from other faces (self-heating)

# Notes
- This is a convenience function that calls all individual flux update functions
- Automatically respects SELF_SHADOWING and SELF_HEATING flags
"""
function update_flux_all!(stpm::SingleAsteroidTPM, r☉::StaticVector{3})
    update_flux_sun!(stpm, r☉)
    update_flux_scat_single!(stpm)
    update_flux_rad_single!(stpm)
end

"""
    update_flux_all!(btpm::BinaryAsteroidTPM, r☉₁::StaticVector{3}, r₁₂::StaticVector{3}, R₁₂::StaticMatrix{3,3})

Update all energy fluxes (solar, scattered, thermal radiation) to the surface for a binary asteroid.
This is a convenience function that computes necessary coordinate transformations and
calls individual flux update functions.

# Arguments
- `btpm::BinaryAsteroidTPM` : Thermophysical model for a binary asteroid
- `r☉₁::StaticVector{3}`    : Sun's position in the primary's body-fixed frame (NOT normalized) [m]
- `r₁₂::StaticVector{3}`    : Position vector of secondary's center in primary's frame [m]
- `R₁₂::StaticMatrix{3,3}`  : Rotation matrix from primary to secondary frame

# Algorithm
1. Computes all necessary coordinate transformations
2. Updates solar flux considering eclipse (mutual shadowing)
3. Updates scattered light flux (self-heating)
4. Updates thermal radiation flux (self-heating)
5. Applies mutual heating between components

# Notes
- This function internally handles all coordinate transformations
- Automatically respects SELF_SHADOWING, SELF_HEATING, MUTUAL_SHADOWING, and MUTUAL_HEATING flags
"""
function update_flux_all!(btpm::BinaryAsteroidTPM, r☉₁::StaticVector{3}, r₁₂::StaticVector{3}, R₁₂::StaticMatrix{3,3})
    # Pre-compute all coordinate transformations
    r☉₂ = R₁₂ * (r☉₁ - r₁₂)  # Sun's position in the secondary's frame
    R₂₁ = R₁₂'               # Rotation matrix from secondary to primary
    r₂₁ = -R₁₂ * r₁₂         # Primary's position in the secondary's frame
    
    # Update all fluxes
    update_flux_sun!(btpm, r☉₁, r☉₂, r₁₂, r₂₁, R₁₂, R₂₁)
    update_flux_scat_single!(btpm)
    update_flux_rad_single!(btpm)
    mutual_heating!(btpm, r₁₂, R₂₁)
end

# ╔═══════════════════════════════════════════════════════════════════╗
# ║                  Energy flux: Sunlight                            ║
# ╚═══════════════════════════════════════════════════════════════════╝


"""
    update_flux_sun!(stpm::SingleAsteroidTPM, r☉::StaticVector{3})

Update the direct solar irradiation flux on every face of the asteroid.

# Arguments
- `stpm::SingleAsteroidTPM` : Thermophysical model for a single asteroid
- `r☉::StaticVector{3}`     : Position vector from asteroid to Sun in body-fixed frame (NOT normalized) [m]

# Algorithm
For each face, the solar flux is calculated as:
1. Solar flux at asteroid's location: F☉ = SOLAR_CONST / distance²
2. Normalize sun direction: r̂☉ = r☉ / |r☉|
3. Face flux: F_sun = F☉ × max(0, n̂ · r̂☉)

where n̂ is the face normal. If `SELF_SHADOWING` is enabled, the function also
checks whether each face is shadowed by other parts of the asteroid.

# Notes
- The input vector `r☉` must not be normalized (used for distance calculation)
- Faces with negative dot product (facing away from Sun) receive zero flux
- Shadowed faces (when `SELF_SHADOWING = true`) also receive zero flux
"""
function update_flux_sun!(stpm::SingleAsteroidTPM, r☉::StaticVector{3})
    # Calculate solar flux and direction
    r̂☉ = normalize(r☉)
    F☉ = SOLAR_CONST / (norm(r☉) * m2au)^2
    
    # Update illumination states
    if stpm.SELF_SHADOWING
        # Check face_visibility_graph availability for self-shadowing
        if isnothing(stpm.shape.face_visibility_graph)
            error(
                "face_visibility_graph must be built when SELF_SHADOWING is enabled. " *
                "Use `build_face_visibility_graph!(shape)` or load shape with `with_face_visibility=true`."
            )
        end
        update_illumination!(stpm.illuminated_faces, stpm.shape, r̂☉; with_self_shadowing=true)
    else
        update_illumination!(stpm.illuminated_faces, stpm.shape, r̂☉; with_self_shadowing=false)
    end

    # Calculate flux for each face
    for i in eachindex(stpm.shape.faces)
        if stpm.illuminated_faces[i]
            n̂ = stpm.shape.face_normals[i]
            stpm.flux_sun[i] = F☉ * (n̂ ⋅ r̂☉)
        else
            stpm.flux_sun[i] = 0.0
        end
    end
end

"""
    update_flux_sun!(
        btpm::BinaryAsteroidTPM,
        r☉₁::StaticVector{3},   r☉₂::StaticVector{3}, 
        r₁₂::StaticVector{3},   r₂₁::StaticVector{3},
        R₁₂::StaticMatrix{3,3}, R₂₁::StaticMatrix{3,3},
    )

Update solar irradiation flux on both components of a binary asteroid system with mutual shadowing.

# Arguments
- `btpm::BinaryAsteroidTPM` : Thermophysical model for a binary asteroid
- `r☉₁::StaticVector{3}`    : Sun's position vector in the primary's body-fixed frame (NOT normalized) [m]
- `r☉₂::StaticVector{3}`    : Sun's position vector in the secondary's body-fixed frame (NOT normalized) [m]
- `r₁₂::StaticVector{3}`    : Position vector of secondary's center in primary's frame [m]
- `r₂₁::StaticVector{3}`    : Position vector of primary's center in secondary's frame [m]
- `R₁₂::StaticMatrix{3,3}`  : Rotation matrix from primary to secondary frame
- `R₂₁::StaticMatrix{3,3}`  : Rotation matrix from secondary to primary frame

# Notes
- All coordinate transformations should be pre-computed by the caller
- Uses the new `apply_eclipse_shadowing!` API from AsteroidShapeModels.jl v0.4.1
- Requires BVH to be built for both shapes (should be done when loading with `with_bvh=true`)
- Combines self-shadowing and mutual shadowing in a single call
"""
function update_flux_sun!(
    btpm::BinaryAsteroidTPM,
    r☉₁::StaticVector{3},   r☉₂::StaticVector{3}, 
    r₁₂::StaticVector{3},   r₂₁::StaticVector{3},
    R₁₂::StaticMatrix{3,3}, R₂₁::StaticMatrix{3,3},
)
    # First, update illumination for both components considering self-shadowing
    update_flux_sun!(btpm.pri, r☉₁)
    update_flux_sun!(btpm.sec, r☉₂)
    
    # Only apply mutual shadowing if enabled
    if btpm.MUTUAL_SHADOWING
        # Check BVH availability
        if isnothing(btpm.pri.shape.bvh) || isnothing(btpm.sec.shape.bvh)
            error(
                "BVH must be built for both shapes when MUTUAL_SHADOWING is enabled. " *
                "Use `build_bvh!(shape)` or load shapes with `with_bvh=true`."
            )
        end

        shape1 = btpm.pri.shape
        shape2 = btpm.sec.shape
        illuminated_faces1 = btpm.pri.illuminated_faces
        illuminated_faces2 = btpm.sec.illuminated_faces

        # Apply eclipse shadowing using the new API from v0.4.1
        # Note: The new API takes the position vector directly instead of translation
        eclipse_status1 = apply_eclipse_shadowing!(illuminated_faces1, shape1, shape2, r☉₁, r₁₂, R₁₂)        
        eclipse_status2 = apply_eclipse_shadowing!(illuminated_faces2, shape2, shape1, r☉₂, r₂₁, R₂₁)
        
        # Update flux_sun based on the updated illumination states
        btpm.pri.flux_sun[.!illuminated_faces1] .= 0.0
        btpm.sec.flux_sun[.!illuminated_faces2] .= 0.0
    end
end

# ╔═══════════════════════════════════════════════════════════════════╗
# ║                 Energy flux: Scattering                           ║
# ╚═══════════════════════════════════════════════════════════════════╝

"""
    update_flux_scat_single!(stpm::SingleAsteroidTPM)

Update flux of scattered sunlight, only considering single scattering.

# Arguments
- `stpm` : Thermophysical model for a single asteroid
"""
function update_flux_scat_single!(stpm::SingleAsteroidTPM)
    stpm.SELF_HEATING == false && return
    isnothing(stpm.shape.face_visibility_graph) && return

    for i_face in eachindex(stpm.shape.faces)
        stpm.flux_scat[i_face] = 0.
        
        # Face properties visible from `i_face`: Face indices and view factors
        visible_indices = get_visible_face_indices(stpm.shape.face_visibility_graph, i_face)
        view_factors = get_view_factors(stpm.shape.face_visibility_graph, i_face)
        
        for (j, fᵢⱼ) in zip(visible_indices, view_factors)
            R_vis = stpm.thermo_params.reflectance_vis[j]

            stpm.flux_scat[i_face] += fᵢⱼ * R_vis * stpm.flux_sun[j]
        end
    end
end

"""
    update_flux_scat_single!(btpm::BinaryAsteroidTPM)

Update flux of scattered sunlight, only considering single scattering.

# Arguments
- `btpm` : Thermophysical model for a binary asteroid
"""
function update_flux_scat_single!(btpm::BinaryAsteroidTPM)
    update_flux_scat_single!(btpm.pri)
    update_flux_scat_single!(btpm.sec)
end

##= TODO: Implement update_flux_scat_mult! =##

# """
#     update_flux_scat_mult!(shape, params::AbstractThermoParams)
#     update_flux_scat_mult!(shape, R_vis)

# Update flux of scattered sunlight, considering multiple scattering.
# """

# ╔═══════════════════════════════════════════════════════════════════╗
# ║                Energy flux: Thermal radiation                     ║
# ╚═══════════════════════════════════════════════════════════════════╝

"""
    update_flux_rad_single!(stpm::SingleAsteroidTPM)

Update flux of absorption of thermal radiation from surrounding surface.
Single radiation-absorption is only considered, assuming albedo is close to zero at thermal infrared wavelength.

# Arguments
- `stpm` : Thermophysical model for a single asteroid
"""
function update_flux_rad_single!(stpm::SingleAsteroidTPM)
    stpm.SELF_HEATING == false && return
    isnothing(stpm.shape.face_visibility_graph) && return

    for i in eachindex(stpm.shape.faces)
        stpm.flux_rad[i] = 0.
        
        # Face properties visible from `i_face`: Face indices and view factors
        visible_indices = get_visible_face_indices(stpm.shape.face_visibility_graph, i)
        view_factors = get_view_factors(stpm.shape.face_visibility_graph, i)
        
        for (j, fᵢⱼ) in zip(visible_indices, view_factors)
            ε    = stpm.thermo_params.emissivity[j]
            R_ir = stpm.thermo_params.reflectance_ir[j]
            Tⱼ   = stpm.temperature[begin, j]
            
            stpm.flux_rad[i] += ε * σ_SB * (1 - R_ir) * fᵢⱼ * Tⱼ^4
        end
    end
end

"""
    update_flux_rad_single!(btpm::BinaryAsteroidTPM)

Update flux of absorption of thermal radiation from surrounding surface.
Single radiation-absorption is only considered, assuming albedo is close to zero at thermal infrared wavelength.

# Arguments
- `btpm` : Thermophysical model for a binary asteroid
"""
function update_flux_rad_single!(btpm::BinaryAsteroidTPM)
    update_flux_rad_single!(btpm.pri)
    update_flux_rad_single!(btpm.sec)
end

# ╔═══════════════════════════════════════════════════════════════════╗
# ║                Mutual heating of binary asteroid                  ║
# ╚═══════════════════════════════════════════════════════════════════╝


"""
    mutual_heating!(btpm::BinaryAsteroidTPM, r₁₂, R₂₁)

Calculate the mutual heating between the primary and secondary asteroids.

# Arguments
- `btpm::BinaryAsteroidTPM` : Thermophysical model for a binary asteroid
- `r₁₂::StaticVector{3}`    : Position vector of secondary's center in primary's frame [m]
- `R₂₁::StaticMatrix{3,3}`  : Rotation matrix from secondary to primary frame

# TODO
- Need to consider local horizon?
"""
function mutual_heating!(btpm::BinaryAsteroidTPM, r₁₂, R₂₁)
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
            c₂ = R₂₁ * c₂ + r₁₂
            n̂₂ = R₂₁ * n̂₂

            f₁₂, d₁₂, d̂₁₂ = view_factor(c₁, c₂, n̂₁, n̂₂, a₂)  # View factor from △A₁B₁C₁ to △A₂B₂C₂
            f₂₁, d₂₁, d̂₂₁ = view_factor(c₂, c₁, n̂₂, n̂₁, a₁)  # View factor from △A₂B₂C₂ to △A₁B₁C₁

            ## if △A₁B₁C₁ and △A₂B₂C₂ are facing each other
            if d̂₁₂ ⋅ n̂₁ > 0 && d̂₁₂ ⋅ n̂₂ < 0
                T₁ = btpm.pri.temperature[begin, i]
                T₂ = btpm.sec.temperature[begin, j]

                ε₁     = thermo_params1.emissivity[i]
                ε₂     = thermo_params2.emissivity[j]
                R_vis₁ = thermo_params1.reflectance_vis[i]
                R_vis₂ = thermo_params2.reflectance_vis[j]
                R_ir₁  = thermo_params1.reflectance_ir[i]
                R_ir₂  = thermo_params2.reflectance_ir[j]

                ## Mutual heating by scattered light
                btpm.pri.flux_scat[i] += f₁₂ * R_vis₂ * btpm.sec.flux_sun[j]
                btpm.sec.flux_scat[j] += f₂₁ * R_vis₁ * btpm.pri.flux_sun[i]

                ## Mutual heating by thermal radiation
                btpm.pri.flux_rad[i] += ε₂ * σ_SB * (1 - R_ir₂) * f₁₂ * T₂^4
                btpm.sec.flux_rad[j] += ε₁ * σ_SB * (1 - R_ir₁) * f₂₁ * T₁^4
            end
        end
    end
end
