# Migration Guide

This page summarizes breaking changes between versions and how to update your code.

---

## v0.3.0

The API changes of v0.2.0 and v0.3.0 are listed in the [changelog](https://github.com/Astroshaper/AsteroidThermoPhysicalModels.jl/blob/main/CHANGELOG.md).

### `ThermoParams` / `GridParams`: material and grid parameters are separate types

```julia
# before
thermo_params = ThermoParams(k, ρ, Cₚ, R_vis, R_ir, ε, z_max, Δz, n_depth)
problem = SingleAsteroidThermoPhysicalProblem(shape, thermo_params; ...)

# after
thermo_params = ThermoParams(; conductivity=k, density=ρ, heat_capacity=Cₚ,
    reflectance_vis=R_vis, reflectance_ir=R_ir, emissivity=ε)
grid_params = GridParams(; z_max, n_depth)   # Δz is auto-computed
problem = SingleAsteroidThermoPhysicalProblem(shape, thermo_params, grid_params; ...)
```

`grid_params` is a required third positional argument of `BinaryAsteroidThermoPhysicalProblem`
as well; a single non-tuple instance is shared by both bodies. Field access changes from
`thermo_params.thermal_conductivity` to `thermo_params.conductivity`, and from
`thermo_params.n_depth` / `.Δz` to `grid_params.n_depth` / `.Δz`. `broadcast_thermo_params!`
is gone: length-1 vectors are expanded to the number of faces when the problem is built.

### `SingleAsteroidOutputSpec`: face-specific outputs are selected by face lists

`subsurface_face_ids` is now a keyword, and the `save_subsurface_temperature` /
`save_roughness_surface_temperature` flags are gone: listing faces switches the output on, an
empty list (the default) leaves it off. Only the surface temperature is saved by default.

```julia
# before
output = SingleAsteroidOutputSpec(output_times, subsurface_face_ids)
output = SingleAsteroidOutputSpec(output_times, Int[]; save_subsurface_temperature=false, save_forces=true)
output = SingleAsteroidOutputSpec(output_times, subsurface_face_ids;
    roughness_face_ids=[1, 7], save_roughness_surface_temperature=true)

# after
output = SingleAsteroidOutputSpec(output_times; subsurface_face_ids)
output = SingleAsteroidOutputSpec(output_times; save_forces=true)
output = SingleAsteroidOutputSpec(output_times; subsurface_face_ids, roughness_face_ids=[1, 7])
```

The `BinaryAsteroidOutputSpec` convenience constructor changes in the same way:
`BinaryAsteroidOutputSpec(output_times_primary, output_times_secondary; subsurface_face_ids_primary, subsurface_face_ids_secondary, ...)`.

### Requires AsteroidShapeModels.jl v0.6

v0.3.0 raises the AsteroidShapeModels.jl compat to `"0.6"`. Surface roughness is carried by
`ShapeModel` itself (`shape.roughness`, built with `add_roughness_models!`); load shapes with
plain `load_shape_obj` / `load_shape_grid`. If your scripts used AsteroidShapeModels.jl v0.5
APIs directly (e.g., `HierarchicalShapeModel`, `as_hierarchical=true`, `has_roughness_model`),
see the [AsteroidShapeModels.jl migration guide](https://astroshaper.github.io/AsteroidShapeModels.jl/stable/guides/migration/).

### Exceptions: input errors raise `ArgumentError`

Errors that you fix by changing the inputs now raise `ArgumentError` instead of a plain
`ErrorException`, as in AsteroidShapeModels.jl v0.6: a missing face visibility graph while
`with_self_shadowing` / `with_self_heating` is enabled, a missing BVH while
`with_mutual_shadowing` is enabled, and an explicit Euler time step at or beyond the stability
limit (λ = αΔt/Δz² ≥ 0.5). Update any `catch` or `@test_throws` that named `ErrorException`.

### Results that change: net thermal force on non-spherical shapes

Up to v0.2.1 the net thermal force was accumulated as ``\sum_i (\hat{\mathbf{r}}_i \cdot \mathbf{F}_i)\,\hat{\mathbf{r}}_i``,
each facet force projected onto the direction of its centre from the origin, instead of
``\sum_i \mathbf{F}_i``. The projection has no physical basis — where a force acts does not enter the
motion of the centre of mass — and dropped the tangential part of every facet force. Corrected in
[#232](https://github.com/Astroshaper/AsteroidThermoPhysicalModels.jl/pull/232).

- **Unaffected**: spheres centred at the origin and symmetric polyhedra whose facet centres lie
  along their normals, where the projection is the identity or its spurious part cancels.
- **Affected**: every irregular shape. As an order of magnitude, on the 49k-facet Ryugu model
  the rotation-averaged net force was about 7 % too small in magnitude and 6° off in direction.

`face_forces` and the torque are unchanged. If you have published or archived net forces
(`forces` in the solution, `thermal_net_forces.csv`) computed with v0.2.1 or earlier on a
non-spherical shape, recompute them with v0.3.0. See *Net force and torque* in the
[physical model](physical_model.md) page.

---

## v0.1.1

No breaking changes. No migration required from v0.1.0.

**Note:** `AsteroidShapeModels.jl` v0.5.x is now supported. If you upgrade to v0.5.x,
the geometry functions `crater_curvature_radius` and `concave_spherical_segment` are no
longer re-exported from this package; use `AsteroidShapeModels` directly instead.
These functions were not part of the public API of this package, so most users are unaffected.

---

## v0.1.0 (from v0.0.7)

### 1. Visibility API changes (via AsteroidShapeModels.jl)

```julia
# Old (v0.0.7)
visible_faces = shape.visiblefacets[face_id]

# New (v0.1.0)
visible_faces = get_visible_face_indices(shape.face_visibility_graph, face_id)
view_factors  = get_view_factors(shape.face_visibility_graph, face_id)
```

### 2. Shape loading

BVH and `face_visibility_graph` are now built automatically when needed by TPM constructors.

```julia
# Both are now equivalent; no need to pre-build manually
shape = load_shape_obj("shape.obj"; scale=1000)

# Optional: pre-build for better performance on large models
shape = load_shape_obj("shape.obj"; scale=1000, with_face_visibility=true, with_bvh=true)
```

### 3. Binary asteroid flux updates

The multiple-call pattern was replaced by a unified API:

```julia
# Old (v0.0.7) — multiple calls required
update_flux_sun!(btpm, r☉₁, r₁₂, R₁₂)
update_flux_scat_single!(btpm)
update_flux_rad_single!(btpm)
mutual_heating!(btpm, r₁₂, R₂₁)

# New (v0.1.0) — unified single call
update_flux_all!(btpm, r☉₁, r₁₂, R₁₂)
```

### 4. Minimum Julia version

Raised from 1.6 to 1.10.
