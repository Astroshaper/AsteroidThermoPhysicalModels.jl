# Roadmap

This document outlines the development roadmap for `AsteroidThermoPhysicalModels.jl`.

## v0.1.0 - Moved geometry processing to `AsteroidShapeModels.jl` (Released 2025-07-13)

The v0.1.0 release marks a significant milestone with stabilized core APIs, critical bug fixes, and improved performance. This release includes full support for `AsteroidShapeModels.jl` v0.4.1 with its eclipse shadowing bug fixes and new unified flux API.

### Completed Tasks ✅

- [x] **Update to AsteroidShapeModels.jl v0.4.1**
  - [x] Updated dependency to v0.4.1 (critical eclipse shadowing bug fixes)
  - [x] Added BVH support to all shape loading
  - [x] Fixed isilluminated API compatibility
  - [x] Updated minimum Julia version to 1.10

- [x] **Refactor illumination API using v0.4.1 new functions**
  - [x] Use new API for self-shadowing with `update_illumination!`
  - [x] Use new API for mutual-shadowing with `apply_eclipse_shadowing!`
  - [x] Implement unified flux update API (`update_flux_all!`)
  - [x] Fix coordinate transformation bug in binary systems

- [x] **Performance improvements**
  - [x] Shadow calculation: ~27x faster (0.379 s/call → 0.014 s/call)
  - [x] Overall component benchmarks ~20x faster for Ryugu
  - [x] Updated and fixed benchmark compatibility

- [x] **Documentation updates**
  - [x] Updated examples to use new type aliases
  - [x] Added comprehensive migration guide in `CHANGELOG.md`
  - [x] Updated benchmark documentation

### Breaking Changes in from v0.0.7 to v0.1.0

1. **Visibility API changes** (via AsteroidShapeModels.jl v0.3.0 → v0.4.1)
- Direct access to `shape.visiblefacets[i]` replaced with function calls
- See migration guide in `CHANGELOG.md` for details

2. **Minimum Julia version**: 1.6 → 1.10

3. **New unified flux API**
- Binary asteroid flux updates now use `update_flux_all!`
- Coordinate transformations centralized

## v0.1.1 - Performance and Test Improvements (Released: 2026-06-08)

- [x] **Update to `AsteroidShapeModels.jl` v0.4.2**
  - Illumination evaluation was optimized in v0.4.2, which will lead to ~2.5x faster shadow calculations.

- [x] **Update to `AsteroidShapeModels.jl` v0.5.x**
  - v0.5.x adds `HierarchicalShapeModel` for surface roughness support and `create_shape_crater`
  - No breaking changes affect this package
  - Updated compat to `"0.4.2, 0.5"` to support both v0.4.x and v0.5.x (v0.4.x support will be dropped in v0.2.0 when `HierarchicalShapeModel` becomes a hard requirement)
  - Removed duplicate geometry functions from `src/roughness.jl` that are now provided by `AsteroidShapeModels.jl`

- [x] **Fix HERA SPICE kernel download URLs**
  - Updated ESA BitBucket URLs to ESAC FTP server (authentication now required on BitBucket)

- [x] **Code organization**
  - [x] Split `src/TPM.jl` into focused files: `solver_types.jl`, `tpm_types.jl`, `tpm_result.jl`, `tpm_run.jl`
  - [x] Remove geometric operation tests that duplicate `AsteroidShapeModels.jl` tests

---
**↓ Planned Releases ↓**
---

## v0.2.0 - API Redesign (Target: 2026)

Redesign the API around a Problem-Solver pattern inspired by `DifferentialEquations.jl`. Separating problem definition from simulation state provides the clean internal architecture needed to implement surface roughness support in v0.3.0.

- [x] **Problem-Solver Architecture**
  - [x] Separate problem definition (shape, thermal parameters, flags) from simulation state (temperatures, fluxes)
  - [x] Implement problem types for single and binary asteroids
  - [x] Adopt `DifferentialEquations.jl`-style workflow: define problem → specify solver → solve

- [x] **Input/Output System Overhaul**
  - [x] Formal ephemerides types with validated fields (replaces informal `NamedTuple`)
  - [x] Optional force and torque output in the inertial frame; requires orientation data in the ephemerides
  - [x] Skipping force/torque computation when orientation data is not provided
  - [x] `OutputSpec` type with Bool flags to independently control each output item (`save_surface_temperature`, `save_subsurface_temperature`, `save_face_forces`, `save_forces`, `save_torques`); replaces individual keyword arguments in `solve`

- [x] **API Cleanup**
  - [x] Remove `subsolar_temperature(r☉, params)` overload; use the explicit scalar form `subsolar_temperature(r☉, R_vis, ε)` instead

## v0.2.1 - Configuration File Support (Target: 2026)

- [ ] **TOML configuration file support** for running parameter surveys without writing Julia code; external data files (shape models, pre-computed ephemerides) are referenced by path within the config
- [ ] **CLI runner**: run simulations from the command line via a configuration file

## v0.2.2 - Solver Quality and Polish (Target: 2026)

- [ ] **Heat conduction solver improvements**
  - [ ] Validate numerical methods against analytical solutions
  - [ ] Optimize implicit solver matrix operations

- [ ] **Binary simulation improvements**
  - [ ] Report eclipse status for each body in binary simulations
  - [ ] Add approximation methods for `mutual_heating!` to reduce computation time

- [ ] **Code quality**
  - [ ] Refactor long heat conduction solver functions
  - [ ] Expand unit test coverage for public functions

## v0.3.0 - Surface Roughness Support (Target: 2026)

Introduce thermophysical modeling of surface roughness using `HierarchicalShapeModel` from `AsteroidShapeModels.jl`, built on the clean Problem-Solver architecture established in v0.2.0. Each global face can optionally carry a roughness model, and an independent mini-TPM is run on its sub-faces (Full Sub-facet TPM). This is the most physically accurate approach and provides a basis for validating simpler approximations in the future.

- [ ] **Roughness-aware problem type**: extend the problem type to accept `HierarchicalShapeModel` and hold independent sub-face state (illumination, flux, temperature, thermal force) for each face

- [ ] **Sub-face flux and temperature calculations**: compute solar flux, self-heating, and 1D heat conduction on sub-faces in their local coordinate frames

- [ ] **Global aggregation**: transform sub-face thermal forces to the global frame and accumulate into body-level force and torque

## v0.4.0 - Performance Optimizations (Target: 2026)

- [ ] **Computational Enhancements**
  - [ ] Multi-threading support
  - [ ] GPU acceleration
  - [ ] Memory optimization

- [ ] **Periodic Simulation Mode**
  - [ ] Support for constant heliocentric distance simulations
  - [ ] Illumination state caching for full rotation
  - [ ] Temperature convergence criteria

## v0.5.0 - Extended Features (Target: 2027)

- [ ] **Extended Physics**
  - [ ] Coupled binary asteroid dynamics
  - [ ] Coupled ejecta/dust/spacecraft dynamics
  - [ ] Sublimation modeling for cometary bodies
  - [ ] Temperature-dependent material properties

## v1.0.0 - Stable Release (Target: 2027)

- [ ] **API stability**
  - [ ] Finalized and stable public API
  - [ ] Comprehensive API documentation
  - [ ] Backwards compatibility guarantees

- [ ] **Complete feature set**
  - [ ] All planned physics models implemented
  - [ ] Validated against observational data

- [ ] **Documentation and tutorials**
  - [ ] Complete user guide
  - [ ] Advanced tutorials for research applications
  - [ ] Performance tuning guide

- [ ] **Quality assurance**
  - [ ] >90% test coverage
  - [ ] Automated performance regression testing
  - [ ] Continuous benchmarking
