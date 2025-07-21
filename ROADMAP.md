# Roadmap

This document outlines the development roadmap for `AsteroidThermoPhysicalModels.jl`.

## Release Process

### For v0.1.1:

1. Update `Project.toml` version to `0.1.1`
2. Update `CHANGELOG.md` with release date
3. Ensure all tests pass
4. Create git tag: `git tag -a v0.1.1 -m "Release v0.1.1"`
5. Push tag: `git push origin v0.1.1`
6. Create GitHub release with detailed notes
7. Register in Julia General Registry
   - Create an issue in the repository with:
     - `@JuliaRegistrator register` command
     - Release notes summarizing key changes
     - If there are breaking changes, clearly document them with "Breaking changes" (not "Breaking Changes") and include migration guide

### Post-release:

1. Update `Project.toml` version to `0.2.0-DEV`
2. Start new "Unreleased" section in `CHANGELOG.md`
3. Update this `ROADMAP.md` for next version

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

---
**↓ Planned Releases ↓**
---

## v0.1.1 - Performance and Test Improvements (Target: August 2025)

- [x] **Update to `AsteroidShapeModels.jl` v0.4.2**
  - Illumination evaluation was optimized in v0.4.2, which will lead to ~2.5x faster shadow calculations.

- [ ] **Code organization and refactoring**
  - [ ] Split long files into smaller, more manageable modules
  - [ ] Refactor long functions for better maintainability (e.g., `implicit_euler!`, `crank_nicolson!`)

- [ ] **mutual_heating! optimization**
  - [ ] Refactor `mutual_heating!` function for clarity
  - [ ] Add approximation methods for faster computation

- [ ] **Enhanced eclipse information**
  - [ ] Modify binary `update_flux_sun!` to return eclipse status
  - [ ] Return `(eclipse_status1, eclipse_status2)` for primary and secondary

- [ ] **Test coverage expansion**
  - [ ] Add comprehensive unit tests for all public functions
  - [ ] Create integration tests for complex workflows
  - [ ] Test edge cases and error conditions
  - [ ] Remove geometric operation tests that duplicate `AsteroidShapeModels.jl` tests

## v0.2.0 - Planned Features (Target: September 2025)

- [ ] **Heat Conduction Solver Enhancements**
  - [ ] Validate numerical methods against analytical solutions
  - [ ] Add accuracy tests for different solvers and boundary conditions
  - [ ] Optimize implicit solver matrix operations
  - [ ] Implement periodic heating benchmarks

- [ ] **Surface roughness integration**
  - [ ] Integrate surface roughness models from `AsteroidShapeModels.jl`
  - [ ] Adapt thermal calculations for sub-facet scale effects
  - [ ] Validate against observations with rough surfaces

## v0.3.0 - API Redesign (Target: October 2025)

- [ ] **Problem-Solver Architecture**
  - [ ] Implement `AbstractThermoPhysicalProblem` interface
  - [ ] Create problem types for single and binary asteroids
  - [ ] Adopt `DifferentialEquations.jl`-style workflow: define problem → specify solver → solve

- [ ] **Input/Output System Overhaul**
  - [ ] Structured input data organization (geometry, thermal, computation parameters)
  - [ ] Configuration file support (TOML/YAML) for parameter surveys
  - [ ] Flexible output system for user-specified quantities
  - [ ] Coordinate transformation support for forces and torques

## v0.4.0 - Performance Optimizations (Target: November 2025)

- [ ] **Computational Enhancements**
  - [ ] Multi-threading support
  - [ ] GPU acceleration
  - [ ] Memory optimization

- [ ] **Periodic Simulation Mode**
  - [ ] Support for constant heliocentric distance simulations
  - [ ] Illumination state caching for full rotation
  - [ ] Temperature convergence criteria

## v0.5.0 - Extended Features (Target: December 2025)

- [ ] **Extended Physics**
  - [ ] Coupled binary asteroid dynamics
  - [ ] Coupled ejecta/dust/spacecraft dynamics
  - [ ] Sublimation modeling for cometary bodies
  - [ ] Temperature-dependent material properties

## v1.0.0 - Stable Release (Target: January 2026)

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
