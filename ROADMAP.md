# Roadmap to v0.1.0

This document outlines the tasks and milestones required for the v0.1.0 release of AsteroidThermoPhysicalModels.jl.

## Release Goals

The v0.1.0 release marks a significant milestone with stabilized core APIs and improved performance tracking capabilities. This release includes breaking changes from the AsteroidShapeModels.jl v0.3.0 update.

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

- [x] **Merge pending PRs**
  - [x] #175 - Remove Format suggestions workflow (Merged)

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

1. **Visibility API changes** (via AsteroidShapeModels.jl v0.3.0)
   - `shape.visiblefacets[i]` → `get_visible_face_indices(shape.face_visibility_graph, i)`
   - Direct view factor access → `get_view_factors(shape.face_visibility_graph, i)`
   - Direct direction access → `get_visible_face_directions(shape.face_visibility_graph, i)`

## Migration Guide

Users upgrading from v0.0.7 to v0.1.0 need to update their code if they directly access visibility data:

```julia
# Old (v0.0.7)
visible_faces = shape.visiblefacets[face_id]

# New (v0.1.0)
visible_faces = get_visible_face_indices(shape.face_visibility_graph, face_id)
```

## Timeline

- **Week 1**: Complete critical tasks and PR merges
- **Week 2**: Documentation and testing updates
- **Week 3**: Final review and release

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