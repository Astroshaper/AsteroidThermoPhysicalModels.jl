# Changelog

All notable changes to this project will be documented in this file.

The format is based on [Keep a Changelog](https://keepachangelog.com/en/1.0.0/),
and this project adheres to [Semantic Versioning](https://semver.org/spec/v2.0.0.html).

## [Unreleased]

### Added
- `illuminated_faces` field to `SingleAsteroidThermoPhysicalModel` for batch illumination processing
- Automatic BVH building in `BinaryAsteroidThermoPhysicalModel` constructor when `MUTUAL_SHADOWING` is enabled
- Automatic face_visibility_graph building in `SingleAsteroidThermoPhysicalModel` constructor when `SELF_SHADOWING` is enabled
- Error checking for required data structures (BVH, face_visibility_graph) in illumination functions
- **Unified flux update API** (#180)
  - New `update_flux_all!` function for both single and binary asteroids
  - Centralizes all coordinate transformations for binary systems
  - Simplifies the interface by handling all flux updates in one call

### Changed
- Updated `AsteroidShapeModels.jl` dependency from v0.2.0 to v0.3.0 (#173)
  - **Breaking**: The visibility API from `AsteroidShapeModels.jl` has changed
  - Direct access to `shape.visiblefacets[i]` is replaced with function calls:
    - `get_visible_face_indices(shape.face_visibility_graph, i)`
    - `get_view_factors(shape.face_visibility_graph, i)`
    - `get_visible_face_directions(shape.face_visibility_graph, i)`
  - Added null checks for `shape.face_visibility_graph` before accessing visibility data
- Refactored illumination calculations to use new batch processing APIs from AsteroidShapeModels.jl v0.4.0
  - `update_flux_sun!` now uses `update_illumination!` for self-shadowing calculations
  - Binary asteroid mutual shadowing now uses `apply_eclipse_shadowing!` API
- Simplified binary asteroid ephemerides structure
  - Removed `sun2` field (secondary's sun position is now computed internally)
  - Removed `S2P` field (inverse transformation is now computed internally)
  - Renamed `sun1` to `sun` for consistency
- Updated to AsteroidShapeModels.jl v0.4.1 (#180)
  - **Critical**: Updated from v0.4.0 to v0.4.1 to get eclipse shadowing bug fixes
  - Migrated to new `apply_eclipse_shadowing!` API that takes position vectors
  - BVH and face_visibility_graph are now built automatically when needed by TPM constructors
  - Shape loading functions support optional pre-building with `with_bvh` and `with_face_visibility`
- **Refactored coordinate transformations** (#180)
  - All coordinate transformations now centralized in `update_flux_all!`
  - `update_flux_sun!` for binary asteroids now accepts pre-computed transformations
  - Improved consistency by renaming `rₛ` to `r₁₂` throughout the codebase
  - Better performance by avoiding duplicate coordinate calculations

### Fixed
- Eclipse shadowing in binary asteroid systems now works correctly with AsteroidShapeModels.jl v0.4.1

### Removed
- Deprecated `update_flux_sun!(btpm, r☉₁, r☉₂)` function (replaced by new API)
- `mutual_shadowing!` function (functionality integrated into `update_flux_sun!`)
- Unused `inverse_transformation` and `transform` functions (coordinate transformations now handled internally)
- Debug scripts used for eclipse shadowing investigation

### Migration Guide from v0.0.7

1. **Visibility API changes** (via AsteroidShapeModels.jl):
   ```julia
   # Old (v0.0.7)
   visible_faces = shape.visiblefacets[face_id]
   
   # New (v0.1.0)
   visible_faces = get_visible_face_indices(shape.face_visibility_graph, face_id)
   view_factors = get_view_factors(shape.face_visibility_graph, face_id)
   ```

2. **Shape loading**:
   ```julia
   # BVH and face_visibility_graph are now built automatically when needed
   shape1 = load_shape_obj("primary.obj"; scale=1000)
   shape2 = load_shape_obj("secondary.obj"; scale=1000)
   
   # Optional: pre-build for better performance
   shape1 = load_shape_obj("primary.obj"; scale=1000, with_face_visibility=true, with_bvh=true)
   ```

3. **Binary asteroid flux updates**:
   ```julia
   # Old (multiple calls)
   update_flux_sun!(btpm, r☉₁, r₁₂, R₁₂)
   update_flux_scat_single!(btpm)
   update_flux_rad_single!(btpm)
   mutual_heating!(btpm, r₁₂, R₂₁)
   
   # New (unified API)
   update_flux_all!(btpm, r☉₁, r₁₂, R₁₂)
   ```

## [0.0.7] - 2025-01-06

### Added
- Implicit solvers for heat conduction equation (#161)
  - `ImplicitEulerSolver`: First-order implicit (backward) Euler method
  - `CrankNicolsonSolver`: Second-order Crank-Nicolson method
  - Both solvers are unconditionally stable
- Functions to calculate thermal radiation flux (#145)
- Store thermal force on every facet in `face_forces` field (#117)
- Basic documentation improvements (#166)

### Changed
- **Breaking**: Renamed abstract and concrete types (#159)
  - Added `AbstractAsteroidThermoPhysicalModel` as the abstract type
  - Added type aliases: `AbstractAsteroidTPM`, `SingleAsteroidTPM`, `BinaryAsteroidTPM`
- **Breaking**: Refactored flux representation in `SingleAsteroidThermoPhysicalModel` (#164)
  - Split `flux` into `flux_sun`, `flux_scat`, and `flux_rad`
  - Renamed `flux_total` to `absorbed_energy_flux` (#167)
- **Breaking**: Unified `UniformThermoParams` and `NonUniformThermoParams` into single `ThermoParams` type (#148)
  - Now supports both uniform and non-uniform thermophysical properties
- **Breaking**: Replaced (`A_TH`, `A_B`) with (`R_ir`, `R_vis`) for reflectance parameters (#143)
- Dimensional depth and time in heat conduction equation (#160)
- Give uniform thermophysical properties to all surface facets automatically (#155)
- Use `NaN` for `E_cons` instead of `missing` (#136)

### Fixed
- Fixed kernel URL in tests (#157)
- Fixed typos in documentation and comments (#165)
- Clipped arguments of inverse trigonometric functions (#126)

### Removed
- Removed dependency on `SPICE.jl` from `src` directory (#158)
- Removed Git dependency in tests (#156)
- Removed `constants.jl` file (#115)

### Dependencies
- Updated MeshIO.jl and GeometryBasics.jl versions (#154)
- Now uses MeshIO.jl for mesh I/O operations (#147)

## [0.0.6] - 2024-03-14

Previous release
