# Changelog

All notable changes to this project will be documented in this file.

The format is based on [Keep a Changelog](https://keepachangelog.com/en/1.0.0/),
and this project adheres to [Semantic Versioning](https://semver.org/spec/v2.0.0.html).

## [Unreleased]

### Added
- `illuminated_faces` field to `SingleAsteroidThermoPhysicalModel` for batch illumination processing
- `inverse_transformation` and `transform` functions for coordinate transformations
- Automatic BVH building in `BinaryAsteroidThermoPhysicalModel` constructor when `MUTUAL_SHADOWING` is enabled
- Automatic face_visibility_graph building in `SingleAsteroidThermoPhysicalModel` constructor when `SELF_SHADOWING` is enabled
- Error checking for required data structures (BVH, face_visibility_graph) in illumination functions

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
  - Removed `S2P` field (inverse transformation is now computed using `inverse_transformation`)
  - Renamed `sun1` to `sun` for consistency
- Updated to AsteroidShapeModels.jl v0.4.0 (#178)
  - Now requires `with_bvh=true` when loading shapes for mutual shadowing
  - Shape loading functions now support BVH building

### Fixed
- TBD

### Removed
- Deprecated `update_flux_sun!(btpm, r☉₁, r☉₂)` function (replaced by new API)
- `mutual_shadowing!` function (functionality integrated into `update_flux_sun!`)

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
