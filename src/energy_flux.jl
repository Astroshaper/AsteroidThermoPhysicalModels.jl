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


# Power absorbed by face `i` of `shape` [W]
function _absorbed_power(state::SingleAsteroidThermoPhysicalState, shape::ShapeModel, i::Integer)
    R_vis  = state.problem.thermo_params.reflectance_vis[i]
    R_ir   = state.problem.thermo_params.reflectance_ir[i]
    F_sun  = state.flux_sun[i]
    F_scat = state.flux_scat[i]
    F_rad  = state.flux_rad[i]

    absorbed_energy_flux(R_vis, R_ir, F_sun, F_scat, F_rad) * shape.face_areas[i]
end

# Power emitted by face `i` of `shape` [W]
function _emitted_power(state::SingleAsteroidThermoPhysicalState, shape::ShapeModel, i::Integer)
    ε = state.problem.thermo_params.emissivity[i]
    T = state.temperature[begin, i]  # Surface temperature

    ε * σ_SB * T^4 * shape.face_areas[i]
end


"""
    integrate_absorbed_power(state::SingleAsteroidThermoPhysicalState) -> Float64

Integrate the absorbed energy flux over all surface facets to obtain total absorbed power [W]:
```
P_abs = Σᵢ F_abs,ᵢ × Aᵢ
```

When the shape carries surface roughness, a face with a roughness model is counted from its
sub-faces and not from the global level: its roughness model is a patch that represents the
face statistically, so the absorbed power of the patch, `Σⱼ F_abs,ⱼ aⱼ` in the units of the
model, is scaled to the area of the face by `Aᵢ / A_proj` (see `update_thermal_force!`).
A face without a roughness model contributes as usual.

# See Also
- `integrate_emitted_power` for the total emitted power
- `absorbed_energy_flux` for the per-facet flux calculation
"""
function integrate_absorbed_power(state::SingleAsteroidThermoPhysicalState)
    shape = state.problem.shape
    isempty(state.roughness_states) && return sum(i -> _absorbed_power(state, shape, i), eachindex(shape.faces))

    P_abs = 0.0
    for (i, k) in enumerate(state.face_roughness_indices)
        if k == 0
            P_abs += _absorbed_power(state, shape, i)
        else
            rs = state.roughness_states[k]
            P_abs += shape.face_areas[i] / projected_area(rs.problem.shape) * integrate_absorbed_power(rs)
        end
    end
    P_abs
end


"""
    integrate_emitted_power(state::SingleAsteroidThermoPhysicalState) -> Float64

Integrate the thermal emission over all surface facets to obtain total emitted power [W]:
```
P_emit = Σᵢ εᵢ × σ × Tᵢ⁴ × Aᵢ
```

In thermal equilibrium, `integrate_emitted_power` ≈ `integrate_absorbed_power`.

When the shape carries surface roughness, a face with a roughness model is counted from its
sub-faces, scaled to the area of the face by `Aᵢ / A_proj`, exactly as in
`integrate_absorbed_power`; the smooth-surface emission of that face is not added.

# See Also
- `integrate_absorbed_power` for the total absorbed power
"""
function integrate_emitted_power(state::SingleAsteroidThermoPhysicalState)
    shape = state.problem.shape
    isempty(state.roughness_states) && return sum(i -> _emitted_power(state, shape, i), eachindex(shape.faces))

    P_emit = 0.0
    for (i, k) in enumerate(state.face_roughness_indices)
        if k == 0
            P_emit += _emitted_power(state, shape, i)
        else
            rs = state.roughness_states[k]
            P_emit += shape.face_areas[i] / projected_area(rs.problem.shape) * integrate_emitted_power(rs)
        end
    end
    P_emit
end


# ╔═══════════════════════════════════════════════════════════════════╗
# ║                     Unified flux update API                       ║
# ╚═══════════════════════════════════════════════════════════════════╝

"""
    update_flux_all!(state::SingleAsteroidThermoPhysicalState, r☉::StaticVector{3})

Update all energy fluxes (solar, scattered, thermal radiation) to the surface for a single asteroid.

# Arguments
- `state` : Thermophysical simulation state for a single asteroid, with or without surface roughness
- `r☉::StaticVector{3}`     : Sun's position in the asteroid-fixed frame (NOT normalized) [m]

# Algorithm
1. Updates direct solar flux on all faces considering self-shadowing
2. Updates scattered sunlight flux from other faces (self-heating)
3. Updates thermal radiation flux from other faces (self-heating)

# Notes
- This is a convenience function that calls all individual flux update functions
- Automatically respects `with_self_shadowing` and `with_self_heating`
- When the shape carries surface roughness, each of the three updates handles the global faces
  first and then the sub-faces of every roughness model. The order of the three calls is
  fixed: the external irradiation of the sub-faces reads the scattered and thermal flux of the
  parent global face, which must therefore be complete before the sub-face update runs.
"""
function update_flux_all!(state::SingleAsteroidThermoPhysicalState, r☉::StaticVector{3})
    update_flux_sun!(state, r☉)
    update_flux_scat_single!(state)
    update_flux_rad_single!(state)
end

"""
    update_flux_all!(state::BinaryAsteroidThermoPhysicalState, r☉₁::StaticVector{3}, r₁₂::StaticVector{3}, R₁₂::StaticMatrix{3,3})

Update all energy fluxes (solar, scattered, thermal radiation) to the surface for a binary asteroid.
This is a convenience function that computes necessary coordinate transformations and
calls individual flux update functions.

# Arguments
- `state::BinaryAsteroidThermoPhysicalState` : Thermophysical simulation state for a binary asteroid
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
function update_flux_all!(state::BinaryAsteroidThermoPhysicalState, r☉₁::StaticVector{3}, r₁₂::StaticVector{3}, R₁₂::StaticMatrix{3,3})
    # Pre-compute all coordinate transformations
    r☉₂ = R₁₂ * (r☉₁ - r₁₂)  # Sun's position in the secondary's frame
    R₂₁ = R₁₂'               # Rotation matrix from secondary to primary
    r₂₁ = -R₁₂ * r₁₂         # Primary's position in the secondary's frame
    
    # Update all fluxes
    update_flux_sun!(state, r☉₁, r☉₂, r₁₂, r₂₁, R₁₂, R₂₁)
    update_flux_scat_single!(state)
    update_flux_rad_single!(state)
    mutual_heating!(state, r₁₂, R₂₁)
end

# ╔═══════════════════════════════════════════════════════════════════╗
# ║                  Energy flux: Sunlight                            ║
# ╚═══════════════════════════════════════════════════════════════════╝


# Raise if `shape` lacks the visibility graph that the modeling flag `flag` requires. The
# problem constructor builds the graph whenever the flag is enabled, so this is reached only
# when the graph was removed from the (mutable) shape afterwards — which must fail loudly
# rather than silently drop the term that depends on it.
function _require_face_visibility_graph(shape::ShapeModel, flag::AbstractString)
    has_face_visibility_graph(shape) || throw(ArgumentError(
        "face_visibility_graph must be built when `$flag` is enabled. " *
        "Use `build_face_visibility_graph!(shape)` or load the shape with `with_face_visibility=true`."
    ))
    return nothing
end


# Shared implementation of the solar flux update for the faces of `shape`. Used for the
# global faces of `state` and, through the sub-face states, for the roughness models.
function _update_flux_sun!(
    state::SingleAsteroidThermoPhysicalState, shape::ShapeModel, r☉::StaticVector{3},
)
    # Calculate solar flux and direction
    r̂☉ = normalize(r☉)
    F☉ = SOLAR_CONST / (norm(r☉) * m2au)^2

    # Update illumination states
    if state.problem.with_self_shadowing
        _require_face_visibility_graph(shape, "with_self_shadowing")
        update_illumination!(state.illuminated_faces, shape, r̂☉; with_self_shadowing=true)
    else
        update_illumination!(state.illuminated_faces, shape, r̂☉; with_self_shadowing=false)
    end

    # Calculate flux for each face
    for i in eachindex(shape.faces)
        if state.illuminated_faces[i]
            n̂ = shape.face_normals[i]
            state.flux_sun[i] = F☉ * (n̂ ⋅ r̂☉)
        else
            state.flux_sun[i] = 0.0
        end
    end
end


"""
    update_flux_sun!(state::SingleAsteroidThermoPhysicalState, r☉::StaticVector{3})

Update the direct solar irradiation flux on every face of the asteroid.

# Arguments
- `state::SingleAsteroidThermoPhysicalState` : Thermophysical simulation state for a single asteroid
- `r☉::StaticVector{3}`     : Position vector from asteroid to Sun in body-fixed frame (NOT normalized) [m]

# Algorithm
For each face, the solar flux is calculated as:
1. Solar flux at asteroid's location: F☉ = SOLAR_CONST / distance²
2. Normalize sun direction: r̂☉ = r☉ / |r☉|
3. Face flux: F_sun = F☉ × max(0, n̂ · r̂☉)

where n̂ is the face normal. If `with_self_shadowing` is enabled, the function also
checks whether each face is shadowed by other parts of the asteroid.

# Notes
- The input vector `r☉` must not be normalized (used for distance calculation)
- Faces with negative dot product (facing away from Sun) receive zero flux
- Shadowed faces (when `with_self_shadowing = true`) also receive zero flux
- When the shape carries surface roughness, the sub-faces of every roughness model are then
  illuminated in the local frame of their parent face (dark when the parent is not illuminated)
"""
function update_flux_sun!(state::SingleAsteroidThermoPhysicalState, r☉::StaticVector{3})
    _update_flux_sun!(state, state.problem.shape, r☉)
    _update_roughness_flux_sun!(state, r☉)
end


# Sub-face solar flux for every face that carries a roughness model; no-op for a smooth
# surface. The global level must be updated first: the gate below reads
# `state.illuminated_faces` — if the parent face is not illuminated (facing away from the Sun,
# or shadowed by the global topography), every one of its sub-faces is dark. A roughness model
# has no terrain beyond its own rim, so without this gate a Sun below the local horizon could
# still light sub-faces tilted towards it. Otherwise the Sun vector is rotated into the local
# frame of the face (a pure rotation, so `|r☉|` and hence the solar flux are preserved) and
# the sub-faces are illuminated like any `ShapeModel`, including self-shadowing within the
# roughness model.
function _update_roughness_flux_sun!(state::SingleAsteroidThermoPhysicalState, r☉::StaticVector{3})
    shape = state.problem.shape
    for (i, k) in enumerate(state.face_roughness_indices)
        k == 0 && continue
        rs = state.roughness_states[k]
        if state.illuminated_faces[i]
            r☉_local = transform_physical_vector_global_to_local(shape, i, r☉)
            update_flux_sun!(rs, r☉_local)
        else
            rs.illuminated_faces .= false
            rs.flux_sun          .= 0.0
        end
    end
end

"""
    update_flux_sun!(
        state::BinaryAsteroidThermoPhysicalState,
        r☉₁::StaticVector{3},   r☉₂::StaticVector{3}, 
        r₁₂::StaticVector{3},   r₂₁::StaticVector{3},
        R₁₂::StaticMatrix{3,3}, R₂₁::StaticMatrix{3,3},
    )

Update solar irradiation flux on both components of a binary asteroid system with mutual shadowing.

# Arguments
- `state::BinaryAsteroidThermoPhysicalState` : Thermophysical simulation state for a binary asteroid
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
    state::BinaryAsteroidThermoPhysicalState,
    r☉₁::StaticVector{3},   r☉₂::StaticVector{3}, 
    r₁₂::StaticVector{3},   r₂₁::StaticVector{3},
    R₁₂::StaticMatrix{3,3}, R₂₁::StaticMatrix{3,3},
)
    # First, update illumination for both components considering self-shadowing
    update_flux_sun!(state.primary, r☉₁)
    update_flux_sun!(state.secondary, r☉₂)
    
    # Only apply mutual shadowing if enabled
    if state.problem.with_mutual_shadowing
        # Check BVH availability
        if !has_bvh(state.primary.problem.shape) || !has_bvh(state.secondary.problem.shape)
            throw(ArgumentError(
                "BVH must be built for both shapes when `with_mutual_shadowing` is enabled. " *
                "Use `build_bvh!(shape)` or load shapes with `with_bvh=true`."
            ))
        end

        shape1 = state.primary.problem.shape
        shape2 = state.secondary.problem.shape
        illuminated_faces1 = state.primary.illuminated_faces
        illuminated_faces2 = state.secondary.illuminated_faces

        # Apply eclipse shadowing using the new API from v0.4.1
        # Note: The new API takes the position vector directly instead of translation
        eclipse_status1 = apply_eclipse_shadowing!(illuminated_faces1, shape1, shape2, r☉₁, r₁₂, R₁₂)        
        eclipse_status2 = apply_eclipse_shadowing!(illuminated_faces2, shape2, shape1, r☉₂, r₂₁, R₂₁)
        
        # Update flux_sun based on the updated illumination states
        state.primary.flux_sun[.!illuminated_faces1] .= 0.0
        state.secondary.flux_sun[.!illuminated_faces2] .= 0.0
    end
end

# ╔═══════════════════════════════════════════════════════════════════╗
# ║                 Energy flux: Scattering                           ║
# ╚═══════════════════════════════════════════════════════════════════╝

# Shared implementation of the single-scattering update for the faces of `shape`. As for
# `_update_flux_sun!`, used for the global faces and, through the sub-face states, for the
# roughness models.
function _update_flux_scat_single!(state::SingleAsteroidThermoPhysicalState, shape::ShapeModel)
    state.problem.with_self_heating == false && return
    _require_face_visibility_graph(shape, "with_self_heating")

    for i_face in eachindex(shape.faces)
        state.flux_scat[i_face] = 0.

        # Face properties visible from `i_face`: Face indices and view factors
        visible_indices = get_visible_face_indices(shape.face_visibility_graph, i_face)
        view_factors = get_view_factors(shape.face_visibility_graph, i_face)

        for (j, fᵢⱼ) in zip(visible_indices, view_factors)
            R_vis = state.problem.thermo_params.reflectance_vis[j]

            state.flux_scat[i_face] += fᵢⱼ * R_vis * state.flux_sun[j]
        end
    end
end


"""
    update_flux_scat_single!(state::SingleAsteroidThermoPhysicalState)

Update flux of scattered sunlight, only considering single scattering.

# Arguments
- `state` : Thermophysical simulation state for a single asteroid
"""
function update_flux_scat_single!(state::SingleAsteroidThermoPhysicalState)
    _update_flux_scat_single!(state, state.problem.shape)

    # Sub-faces of every roughness model (empty loop for a smooth surface): scattering between
    # the sub-faces, plus the sunlight reflected towards them by the other global faces — from
    # the sub-faces of those faces' own roughness models, direction by direction (see
    # `_add_external_scattering!`). Only with self-heating, which is what the external term is.
    for (i, k) in enumerate(state.face_roughness_indices)
        k == 0 && continue
        rs = state.roughness_states[k]
        _update_flux_scat_single!(rs, rs.problem.shape)
        state.problem.with_self_heating && _add_external_scattering!(state, k, i)
    end
end


"""
    _directional_emission(patch::ShapeModel, visible, d̂_local, emission) -> E

Emission of a roughness model `patch` towards the direction `d̂_local` (unit vector in the local
frame of the patch), per unit area of the patch's reference plane projected onto that direction:
```
E(d̂) = Σₙ Vₙ (n̂ₙ ⋅ d̂)⁺ Eₙ aₙ / (A_proj cos θ)
```
where `visible[n]` (`Vₙ`) tells whether sub-face `n` is seen from `d̂`, `emission(n)` (`Eₙ`)
is the quantity emitted by sub-face `n` per unit area — `ε σ T⁴` for thermal emission,
`R_vis F_sun` for reflected sunlight, a radiance for `roughness_radiance` — and `cos θ` is the
z component of `d̂_local`. For a uniform, unshadowed patch this reduces to `E = Eₙ`: the
representative patch radiates like a smooth Lambertian face. Sub-faces facing away from `d̂`
do not contribute.

The caller guarantees `cos θ > 0`.
"""
function _directional_emission(patch::ShapeModel, visible, d̂_local::StaticVector{3}, emission)
    cosθ  = d̂_local[3]
    total = 0.0
    for n in eachindex(patch.faces)
        visible[n] || continue
        cosθ_n = patch.face_normals[n] ⋅ d̂_local
        cosθ_n ≤ 0 && continue
        total += cosθ_n * patch.face_areas[n] * emission(n)
    end
    return total / (projected_area(patch) * cosθ)
end


# Add to the sub-face fluxes `flux_sub` of `patch` the irradiance `F` that the parent face
# receives from the direction `d̂_local` (local frame): a distant source of intensity `F / cos θ`
# lights the sub-faces it can see, each by its own inclination — the same rule as for the
# sunlight. The power received, Σₘ Fₘ aₘ, equals F A_proj up to the shadowing discretisation
# of the sub-faces (seen or not from the ray through their centre): on a height field every
# ray through the reference plane meets a sub-face, the walls shading the floor.
function _distribute_directional!(flux_sub::AbstractVector, patch::ShapeModel, visible, d̂_local::StaticVector{3}, F::Real)
    cosθ = d̂_local[3]
    (F ≤ 0 || cosθ ≤ 0) && return
    intensity = F / cosθ
    for m in eachindex(patch.faces)
        visible[m] || continue
        cosθ_m = patch.face_normals[m] ⋅ d̂_local
        cosθ_m > 0 && (flux_sub[m] += intensity * cosθ_m)
    end
end


"""
    _add_external_radiation!(state, k, i)
    _add_external_scattering!(state, k, i)

Add to the sub-faces of the roughness model on global face `i` (sub-state `k`) the thermal
radiation, respectively the reflected sunlight, that they receive from the other global faces.

For every face `j` visible from `i` (view factor `f_ij`, direction `d̂_ij` from the face
visibility graph):

1. **Emitter**: the emission of `j` towards `i`. If `j` carries a roughness model, it is the
   directional emission of that model (`_directional_emission`) — the hot sunlit wall of a
   crater that faces `i` radiates more towards `i` than a smooth face would (thermal-infrared
   beaming), a shaded wall less. The sub-faces of `j` seen from `i` come from the mask that `j`
   precomputed towards `i` (`RoughnessNeighbours`). If `j` is smooth, it radiates as a Lambertian
   face, `ε σ T_j⁴` or `R_vis F_sun,j`.
2. **Far-field**: the irradiance reaching face `i` is `F_ij = E_j(d̂_ji) f_ij`, the same form as
   the smooth-face term `ε σ T_j⁴ f_ij` — the view factor already carries the geometry of the
   pair, and both faces are taken as small compared to their distance, as the view factor does.
3. **Receiver**: `F_ij` is distributed over the sub-faces of `i` that see `j`
   (`_distribute_directional!`), so the wall facing `j` is heated and the wall behind it is not.

Thermal emission uses the sub-face temperatures of the previous time step and reflected sunlight
uses the direct solar flux only (single scattering, as on the global level), so the result does
not depend on the order in which the roughness models are updated. The global-level fluxes of
face `i` are not touched: they remain the smooth-surface baseline.
"""
function _add_external_radiation!(state::SingleAsteroidThermoPhysicalState, k::Integer, i::Integer)
    _add_external!(state, k, i, :rad)
end

function _add_external_scattering!(state::SingleAsteroidThermoPhysicalState, k::Integer, i::Integer)
    _add_external!(state, k, i, :scat)
end

function _add_external!(state::SingleAsteroidThermoPhysicalState, k::Integer, i::Integer, kind::Symbol)
    shape      = state.problem.shape
    graph      = shape.face_visibility_graph
    rs         = state.roughness_states[k]
    neighbours = state.roughness_neighbours[k]
    patch      = rs.problem.shape
    flux_sub   = kind === :rad ? rs.flux_rad : rs.flux_scat

    visible_indices = get_visible_face_indices(graph, i)
    view_factors    = get_view_factors(graph, i)
    directions      = get_visible_face_directions(graph, i)

    for (p, (j, fᵢⱼ, d̂ᵢⱼ)) in enumerate(zip(visible_indices, view_factors, directions))
        k_j = state.face_roughness_indices[j]

        # Emitter: j towards i
        if k_j == 0
            Eⱼ = kind === :rad ?
                state.problem.thermo_params.emissivity[j] * σ_SB * state.temperature[begin, j]^4 :
                state.problem.thermo_params.reflectance_vis[j] * state.flux_sun[j]
        else
            rs_j    = state.roughness_states[k_j]
            patch_j = rs_j.problem.shape
            visible_sub_faces_j = state.roughness_neighbours[k_j].visible_sub_faces[neighbours.position_in_neighbour[p]]
            d̂ⱼᵢ     = transform_physical_vector_global_to_local(shape, j, -d̂ᵢⱼ)
            tp_j    = rs_j.problem.thermo_params
            Eⱼ = kind === :rad ?
                _directional_emission(patch_j, visible_sub_faces_j, d̂ⱼᵢ, n -> tp_j.emissivity[n] * σ_SB * rs_j.temperature[begin, n]^4) :
                _directional_emission(patch_j, visible_sub_faces_j, d̂ⱼᵢ, n -> tp_j.reflectance_vis[n] * rs_j.flux_sun[n])
        end

        # Far-field: irradiance on face i, then distributed over the sub-faces that see j
        Fᵢⱼ = Eⱼ * fᵢⱼ
        d̂ᵢⱼ_local = transform_physical_vector_global_to_local(shape, i, d̂ᵢⱼ)
        _distribute_directional!(flux_sub, patch, neighbours.visible_sub_faces[p], d̂ᵢⱼ_local, Fᵢⱼ)
    end
end


# Sky view factor of face `j`: the part of its hemisphere not covered by the other faces of
# the same shape, `1 − Σₖ fⱼₖ`. Clamped at zero against round-off in the view factors.
_sky_view_factor(graph::FaceVisibilityGraph, j::Integer) = max(0.0, 1 - sum(get_view_factors(graph, j)))

"""
    update_flux_scat_single!(state::BinaryAsteroidThermoPhysicalState)

Update flux of scattered sunlight, only considering single scattering.

# Arguments
- `state` : Thermophysical simulation state for a binary asteroid
"""
function update_flux_scat_single!(state::BinaryAsteroidThermoPhysicalState)
    update_flux_scat_single!(state.primary)
    update_flux_scat_single!(state.secondary)
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

# Shared implementation of the thermal-radiation update for the faces of `shape`. `flux_rad`
# is the irradiance incident on each face from the others — like `flux_sun` and `flux_scat`,
# it carries no absorptivity: the receiving face's `1 − R_ir` is applied where the flux is
# absorbed (surface boundary condition, `absorbed_energy_flux`) and its `R_ir` where the
# reflected part is emitted (`update_thermal_force!`).
function _update_flux_rad_single!(state::SingleAsteroidThermoPhysicalState, shape::ShapeModel)
    state.problem.with_self_heating == false && return
    _require_face_visibility_graph(shape, "with_self_heating")

    for i in eachindex(shape.faces)
        state.flux_rad[i] = 0.

        # Face properties visible from `i_face`: Face indices and view factors
        visible_indices = get_visible_face_indices(shape.face_visibility_graph, i)
        view_factors = get_view_factors(shape.face_visibility_graph, i)

        for (j, fᵢⱼ) in zip(visible_indices, view_factors)
            εⱼ = state.problem.thermo_params.emissivity[j]
            Tⱼ = state.temperature[begin, j]

            # Irradiance on i from the Lambertian emission of j: Eⱼ Aⱼ Fⱼᵢ / Aᵢ = Eⱼ Fᵢⱼ by
            # reciprocity, with Fᵢⱼ the view factor stored for face i
            state.flux_rad[i] += εⱼ * σ_SB * fᵢⱼ * Tⱼ^4
        end
    end
end


"""
    update_flux_rad_single!(state::SingleAsteroidThermoPhysicalState)

Update the flux of thermal radiation incident on each face from the surrounding surface.

Only the direct emission of the other faces is counted (single bounce); thermal radiation they
reflect is neglected. `flux_rad` is an incident flux, like `flux_sun` and `flux_scat`: the
thermal-infrared reflectance of the receiving face is applied where the flux is absorbed.

# Arguments
- `state` : Thermophysical simulation state for a single asteroid
"""
function update_flux_rad_single!(state::SingleAsteroidThermoPhysicalState)
    _update_flux_rad_single!(state, state.problem.shape)

    # Sub-faces of every roughness model (empty loop for a smooth surface): radiative exchange
    # between the sub-faces from their surface temperatures, plus the thermal radiation the
    # other global faces send towards them — from the sub-faces of those faces' own roughness
    # models, direction by direction (see `_add_external_radiation!`). Only with self-heating.
    for (i, k) in enumerate(state.face_roughness_indices)
        k == 0 && continue
        rs = state.roughness_states[k]
        _update_flux_rad_single!(rs, rs.problem.shape)
        state.problem.with_self_heating && _add_external_radiation!(state, k, i)
    end
end

"""
    update_flux_rad_single!(state::BinaryAsteroidThermoPhysicalState)

Update flux of absorption of thermal radiation from surrounding surface.
Single radiation-absorption is only considered, assuming albedo is close to zero at thermal infrared wavelength.

# Arguments
- `state` : Thermophysical simulation state for a binary asteroid
"""
function update_flux_rad_single!(state::BinaryAsteroidThermoPhysicalState)
    update_flux_rad_single!(state.primary)
    update_flux_rad_single!(state.secondary)
end

# ╔═══════════════════════════════════════════════════════════════════╗
# ║                Mutual heating of binary asteroid                  ║
# ╚═══════════════════════════════════════════════════════════════════╝


"""
    mutual_heating!(state::BinaryAsteroidThermoPhysicalState, r₁₂, R₂₁)

Calculate the mutual heating between the primary and secondary asteroids.

# Arguments
- `state::BinaryAsteroidThermoPhysicalState` : Thermophysical simulation state for a binary asteroid
- `r₁₂::StaticVector{3}`    : Position vector of secondary's center in primary's frame [m]
- `R₂₁::StaticMatrix{3,3}`  : Rotation matrix from secondary to primary frame

# TODO
- Need to consider local horizon?
"""
function mutual_heating!(state::BinaryAsteroidThermoPhysicalState, r₁₂, R₂₁)
    state.problem.with_mutual_heating == false && return

    shape1 = state.primary.problem.shape
    shape2 = state.secondary.problem.shape
    thermo_params1 = state.primary.problem.thermo_params
    thermo_params2 = state.secondary.problem.thermo_params

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
                T₁ = state.primary.temperature[begin, i]
                T₂ = state.secondary.temperature[begin, j]

                ε₁     = thermo_params1.emissivity[i]
                ε₂     = thermo_params2.emissivity[j]
                R_vis₁ = thermo_params1.reflectance_vis[i]
                R_vis₂ = thermo_params2.reflectance_vis[j]

                ## Mutual heating by scattered light
                state.primary.flux_scat[i] += f₁₂ * R_vis₂ * state.secondary.flux_sun[j]
                state.secondary.flux_scat[j] += f₂₁ * R_vis₁ * state.primary.flux_sun[i]

                ## Mutual heating by thermal radiation (incident flux; the receiving face's
                ## thermal-infrared reflectance is applied where the flux is absorbed)
                state.primary.flux_rad[i] += ε₂ * σ_SB * f₁₂ * T₂^4
                state.secondary.flux_rad[j] += ε₁ * σ_SB * f₂₁ * T₁^4
            end
        end
    end
end
