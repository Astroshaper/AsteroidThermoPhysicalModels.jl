# Roadmap to v0.1.0

This document outlines the tasks and milestones required for the v0.1.0 release of AsteroidThermoPhysicalModels.jl.

## Release Goals

The v0.1.0 release marks a significant milestone with stabilized core APIs and improved performance tracking capabilities. This release includes updates to support AsteroidShapeModels.jl v0.4.0 with its latest improvements.

## Required Tasks

### 🔴 Critical (Must Have)

- [ ] **Update version number**
  - [ ] Change `Project.toml` version from `0.0.8-DEV` to `0.1.0-DEV`
  - [ ] Final version bump to `0.1.0` before release

- [ ] **Update CHANGELOG.md**
  - [ ] Move Unreleased items to v0.1.0 section
  - [ ] Add release date
  - [ ] Ensure all breaking changes are clearly documented
  - [ ] Add migration guide for visibility API changes

- [x] **Update to AsteroidShapeModels.jl v0.4.1** (Completed)
  - [x] Updated dependency to v0.4.1 (critical eclipse shadowing bug fixes)
  - [x] Added BVH support to all shape loading
  - [x] Fixed isilluminated API compatibility
  - [x] Updated minimum Julia version to 1.10

- [x] **Refactor illumination API using `AsteroidShapeModels.jl` v0.4.1 new functions**
  - [x] Use new API for self-shadowing: Replace current illumination calculations with `update_illumination!`
  - [x] Use new API for mutual-shadowing: Migrate to new `apply_eclipse_shadowing!` API that takes position vectors
  - [x] Implement unified flux update API (`update_flux_all!`) for cleaner interface

- [x] **Merge pending PRs**
  - [x] #175 - Remove Format suggestions workflow (Merged)
  - [x] #177 - Update benchmark documentation (Merged)
  - [x] #178 - Update to AsteroidShapeModels.jl v0.4.0 (Merged)

### 🟡 Important (Should Have)

- [ ] **Documentation updates**
  - [ ] Update examples to use new visibility API
  - [ ] Add migration guide for users upgrading from v0.0.7
  - [ ] Update API documentation for changed functions

- [ ] **Testing improvements**
  - [ ] Add tests for new visibility API usage
  - [ ] Ensure all tests pass with AsteroidShapeModels v0.3.0

### 🟢 Nice to Have (Could Have)

- [ ] **Performance optimizations**
  - [ ] Review benchmark results and identify optimization opportunities
  - [ ] Consider implementing parallel processing for shadow calculations

- [ ] **Code quality**
  - [ ] Refactor long functions identified in code analysis
    - [ ] `mutual_shadowing!` (127 lines)
    - [ ] `implicit_euler!` (74 lines)
    - [ ] `crank_nicolson!` (80 lines)

## Breaking Changes Summary

### From v0.0.7 to v0.1.0

- Visibility API changes via AsteroidShapeModels.jl v0.3.0 → v0.4.1
- Minimum Julia version: 1.6 → 1.10

See CHANGELOG.md for detailed migration guide.

## Future Versions

### v0.2.0 (Planned)
- [ ] Implement surface roughness models
- [ ] Add multi-threading support
- [ ] Improve test coverage for `energy_flux.jl` and `non_grav.jl`

### v1.0.0 (Long-term)
- [ ] Stable API with all major features implemented
- [ ] Comprehensive documentation and tutorials
- [ ] Performance optimizations completed
- [ ] Full test coverage achieved

## Release Checklist

Before releasing v0.1.0:

- [ ] All tests pass locally
- [ ] CI/CD pipeline is green
- [ ] Documentation builds successfully
- [ ] CHANGELOG.md is updated
- [ ] Version number is correct in Project.toml
- [ ] Tag is created and pushed
- [ ] Release notes are written
- [ ] Package is registered/updated in Julia General Registry