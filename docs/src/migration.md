# Migration Guide

This page summarizes breaking changes between versions and how to update your code.

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
