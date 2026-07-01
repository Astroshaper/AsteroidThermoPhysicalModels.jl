# Changelog

All notable changes to this project will be documented in this file.

The format is based on [Keep a Changelog](https://keepachangelog.com/en/1.0.0/),
and this project adheres to [Semantic Versioning](https://semver.org/spec/v2.0.0.html).

## [Unreleased]

## [0.3.0] - TBD

This release separates thermophysical material properties (`ThermoParams`) from
numerical grid configuration (`GridParams`). The old monolithic `ThermoParams` that
held both material and grid parameters is replaced by two focused types.

### Migration Guide

```julia
# v0.2.x
thermo_params = ThermoParams(k, ρ, Cₚ, R_vis, R_ir, ε, z_max, Δz, n_depth)
problem = SingleAsteroidThermoPhysicalProblem(shape, thermo_params; ...)

# v0.3.0 — material and grid parameters are now separate
thermo_params = ThermoParams(
    conductivity    = k,
    density         = ρ,
    heat_capacity   = Cₚ,
    reflectance_vis = R_vis,
    reflectance_ir  = R_ir,
    emissivity      = ε,
)
grid_params = GridParams(; z_max=z_max, n_depth=n_depth)  # Δz is auto-computed
problem = SingleAsteroidThermoPhysicalProblem(shape, thermo_params, grid_params; ...)
```

For binary asteroid problems:

```julia
# v0.2.x
problem = BinaryAsteroidThermoPhysicalProblem(
    (shape1, shape2),
    (thermo_params1, thermo_params2);
    ...
)

# v0.3.0 — grid_params is now a required third positional argument
problem = BinaryAsteroidThermoPhysicalProblem(
    (shape1, shape2),
    (thermo_params1, thermo_params2),
    (grid_params1, grid_params2);
    ...
)
# Or pass a single shared instance for both bodies:
problem = BinaryAsteroidThermoPhysicalProblem(
    (shape1, shape2),
    thermo_params,
    grid_params;
    ...
)
```

For code that reads fields of `ThermoParams`:

```julia
# v0.2.x
thermo_params.thermal_conductivity  # Vector{Float64}
thermo_params.n_depth               # Int
thermo_params.Δz                    # Float64

# v0.3.0
thermo_params.conductivity  # renamed field
grid_params.n_depth
grid_params.Δz
```

### Added

- **`GridParams`** struct for numerical depth-grid configuration, extracted from the
  old monolithic `ThermoParams`:
  - Fields: `z_max` (depth of lower boundary [m]), `n_depth` (number of depth nodes), `Δz` (node spacing [m])
  - `GridParams(; z_max, n_depth)` keyword constructor: `Δz` is auto-computed as `z_max / (n_depth - 1)`, placing nodes uniformly at `0, Δz, 2Δz, …, z_max`; use the positional constructor `GridParams(z_max, n_depth, Δz)` to specify `Δz` explicitly
  - `GridParams` is now exported
- **`ThermoParams` keyword constructor** `ThermoParams(; conductivity, density, heat_capacity, reflectance_vis, reflectance_ir, emissivity)`: avoids relying on positional argument order
- **`ThermoParams` mixed scalar/vector constructor**: scalar arguments are automatically broadcast to match the length of any vector arguments, removing the need for manual `fill()` calls in non-uniform surface cases

### Changed

- **Breaking**: `ThermoParams` no longer holds grid parameters (`z_max`, `Δz`, `n_depth`); pass a `GridParams` instance separately
- **Breaking**: `ThermoParams.thermal_conductivity` renamed to `ThermoParams.conductivity`
- **Breaking**: Positional 9-argument constructor `ThermoParams(k, ρ, Cₚ, R_vis, R_ir, ε, z_max, Δz, n_depth)` removed; use `ThermoParams(k, ρ, Cₚ, R_vis, R_ir, ε)` (6-arg) or the keyword constructor
- **Breaking**: `SingleAsteroidThermoPhysicalProblem(shape, thermo_params; ...)` → `SingleAsteroidThermoPhysicalProblem(shape, thermo_params, grid_params; ...)`
- **Breaking**: `BinaryAsteroidThermoPhysicalProblem((shape1, shape2), (tp1, tp2); ...)` → `BinaryAsteroidThermoPhysicalProblem((shape1, shape2), (tp1, tp2), (gp1, gp2); ...)`; a single `ThermoParams`/`GridParams` instance (non-tuple) can be shared between both bodies

### Removed

- **Breaking**: `broadcast_thermo_params!` removed; single-body `ThermoParams` (length-1 vectors) are now expanded to `n_face` at `SingleAsteroidThermoPhysicalProblem` construction time via an internal `_expand_thermo_params` call, eliminating the need for mutation

### Internal

- Moved `GridParams` definition to `src/grid_params.jl`

## [0.2.1] - 2026-06-26

No migration required from v0.2.0.

### Added

- `*Ephemerides` constructors now accept `AbstractRange` for `times` (e.g. `range(et_begin, et_end; length=n)`); the range is collected to `Vector{Float64}` internally, eliminating the manual `collect` call
- `*Ephemerides` constructors now auto-convert plain `Vector` inputs to `Vector{SVector{3,Float64}}` for position vectors and `Vector{SMatrix{3,3,Float64,9}}` for rotation matrices; no need to import `StaticArrays` at the call site

### Internal

- Renamed `energy_in` / `energy_out` fields to `absorbed_power` / `emitted_power` in the diagnostics data — aligns the implementation with the documented names; these fields are not exported

## [0.2.0] - 2026-06-24

This release introduces a Problem-Solver API redesign inspired by `DifferentialEquations.jl`.
The old `run_TPM!` function is replaced by `solve(problem, algorithm; kwargs...)`.

### Migration Guide

```julia
# v0.1.x
ephem = (time = times, sun = r_sun)   # NamedTuple with `time` and `sun` fields
stpm  = SingleAsteroidThermoPhysicalModel(shape, thermo_params; ...)
result = run_TPM!(stpm, ephem, times_to_save, face_ID)

# v0.2.0
ephem   = SingleAsteroidEphemerides(times, r_sun)
# ephem = SingleAsteroidEphemerides(times, r_sun, R_body_to_inertial)  # for force/torque output

problem  = SingleAsteroidThermoPhysicalProblem(shape, thermo_params; ...)
output   = SingleAsteroidOutputSpec(output_times, subsurface_face_ids;
    save_surface_temperature    = true,
    save_subsurface_temperature = true,
    save_face_forces            = false,
    save_forces                 = false,  # true requires R_body_to_inertial in ephem
    save_torques                = false,  # true requires R_body_to_inertial in ephem
)
solution = solve(problem, CrankNicolson();
    ephem               = ephem,
    output              = output,
    initial_temperature = 200.0,
)
```

### Added

- **Problem-Solver API** via `CommonSolve.jl`
  - `SingleAsteroidThermoPhysicalProblem(shape, thermo_params; kwargs...)`: encapsulates the physical problem definition for a single asteroid
  - `BinaryAsteroidThermoPhysicalProblem(primary, secondary; kwargs...)`: encapsulates the physical problem for a binary asteroid system
  - `BinaryAsteroidThermoPhysicalProblem(shape, thermo_params; kwargs...)`: convenience constructor accepting tuples `(shape1, shape2)` and `(thermo_params1, thermo_params2)`; single-body kwargs are applied identically to both bodies
  - `solve(problem, algorithm; ephem, output, initial_temperature, ...)`: run a simulation for a single asteroid
  - `solve(problem, algorithm; ephem, output, initial_temperature_primary, initial_temperature_secondary, ...)`: run a simulation for a binary system

- **Output specification types**
  - `SingleAsteroidOutputSpec(output_times, subsurface_face_ids; save_surface_temperature=true, save_subsurface_temperature=true, save_face_forces=false, save_forces=false, save_torques=false)`: specifies which timesteps, face indices, and physical quantities to record; Bool flags independently control each output item
  - `BinaryAsteroidOutputSpec(output_times_primary, subsurface_face_ids_primary, output_times_secondary, subsurface_face_ids_secondary; save_*=...)`: convenience constructor with shared Bool flags for both bodies
  - `BinaryAsteroidOutputSpec(primary, secondary)`: base constructor accepting an independent `SingleAsteroidOutputSpec` per body

- **Algorithm types** (previously flags in `run_TPM!`)
  - `ExplicitEuler()`, `ImplicitEuler()`, `CrankNicolson()`

- **Boundary condition types** (previously flags in `run_TPM!`)
  - `RadiationBoundaryCondition()`, `InsulationBoundaryCondition()`, `IsothermalBoundaryCondition()`

- **Ephemerides types** — inputs for `solve`
  - `SingleAsteroidEphemerides(times, r_sun)`: temperature-only mode; `r_sun` is the sun position in the body-fixed frame
  - `SingleAsteroidEphemerides(times, r_sun, R_body_to_inertial)`: force/torque mode; additionally accepts body-to-inertial rotation matrices
  - `BinaryAsteroidEphemerides(times, r_sun, r_secondary, R_primary_to_secondary)`: temperature-only mode for binary systems
  - `BinaryAsteroidEphemerides(times, r_sun, r_secondary, R_primary_to_secondary, R_primary_to_inertial)`: force/torque mode for binary systems
  - The presence or absence of rotation matrices is encoded in the type parameter `{R}` (`R = Nothing` or `R <: AbstractVector`), enabling zero-overhead dispatch in `solve`

- **Solution types**
  - `SingleAsteroidThermoPhysicalSolution`: holds simulation output for a single asteroid; optional fields (`surface_temperature`, `subsurface_temperature`, `face_forces`, `forces`, `torques`) are `nothing` when the corresponding `save_*` flag is `false`
  - `BinaryAsteroidThermoPhysicalSolution`: holds simulation output for a binary asteroid system; contains `primary` and `secondary` sub-solutions

- `export_solution(dirpath, solution)`: export simulation results to CSV files; replaces `export_TPM_results`; files written depend on `output` flags:
  - `diagnostics.csv`: always (`absorbed_power`, `emitted_power` at all timesteps)
  - `surface_temperature.csv`: when `save_surface_temperature = true`
  - `subsurface_temperature.csv`: when `save_subsurface_temperature = true`
  - `thermal_face_forces.csv`: when `save_face_forces = true`
  - `thermal_net_forces.csv`: when `save_forces = true` or `save_torques = true`
- `ThermoParams` is now publicly exported (previously required `AsteroidThermoPhysicalModels.ThermoParams`)

### Changed

- Binary asteroid output directory names changed from `pri`/`sec` to `primary`/`secondary`

### Removed

- `run_TPM!`: replaced by `solve(problem, algorithm; kwargs...)`
- `export_TPM_results`: replaced by `export_solution`
- `subsolar_temperature(r☉, params::AbstractThermoParams)` overload: used `params.reflectance_vis[begin]` and `params.emissivity[begin]` silently, which is semantically inconsistent for non-uniform `ThermoParams`; use `subsolar_temperature(r☉, R_vis, ε)` directly instead

### Internal

- Introduced `SingleAsteroidThermoPhysicalState` / `BinaryAsteroidThermoPhysicalState` as internal simulation state types (separate from problem definition; not exported)
- `init_temperature!` is no longer exported; use the `initial_temperature` keyword argument of `solve` to set the initial temperature
- Renamed `src/tpm_types.jl` → `src/tpm_state.jl`
- Renamed `src/tpm_result.jl` → `src/tpm_solution.jl`
- Renamed `src/tpm_run.jl` → `src/tpm_init.jl` (file now contains only `init_temperature!`)
- Moved `subsolar_temperature` from `tpm_run.jl` to `thermo_params.jl`
- Dropped `AsteroidShapeModels.jl` v0.4.x compat; v0.5+ is now required

## [0.1.1] - 2026-06-08

No migration required from v0.1.0. There are no breaking changes in this release.

### Fixed
- Updated HERA SPICE kernel download URLs to ESAC FTP server in `test/TPM_Didymos/TPM_Didymos.jl` and `benchmark/benchmarks.jl` (#191)
  - The ESA BitBucket server now requires authentication, which caused CI to fail when downloading kernels

### Changed
- Updated `AsteroidShapeModels.jl` dependency from v0.4.1 to v0.4.2
  - Shadow calculations are now 2.7x faster (1.04s → 0.39s for 72 time steps and 49k Ryugu shape)
  - Overall simulation performance improved by 13.6% for single asteroids (Ryugu: 103s → 89s for 20 rotations, 1440 time steps)
  - Face maximum elevation optimization automatically applied when using `with_face_visibility=true`
  - No code changes required due to backward compatibility
- Extended `AsteroidShapeModels.jl` compat to `"0.4.2, 0.5"` to support v0.5.x (#194)
  - v0.5.x adds `HierarchicalShapeModel` and `create_shape_crater`, which are prerequisites for surface roughness support planned in v0.2.0
  - Compatibility with v0.4.x is maintained; v0.4.x support will be dropped in v0.2.0

### Removed
- Removed `crater_curvature_radius` and `concave_spherical_segment` from `src/roughness.jl`; these functions are now provided by `AsteroidShapeModels.jl` v0.5.0+ (#192)
- Removed `test/find_visiblefacets.jl`; the test duplicated coverage already provided by `AsteroidShapeModels.jl` (#195)

### Internal
- Split `src/TPM.jl` (892 lines) into focused files: `solver_types.jl`, `tpm_types.jl`, `tpm_result.jl`, `tpm_run.jl` (#195)

## [0.1.0] - 2025-07-13

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
- Coordinate transformation bug in binary systems where `r₂₁` was incorrectly calculated (#183)

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
