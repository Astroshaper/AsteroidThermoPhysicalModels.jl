# AsteroidThermoPhysicalModels.jl

[![Stable](https://img.shields.io/badge/docs-stable-blue.svg)](https://Astroshaper.github.io/AsteroidThermoPhysicalModels.jl/stable)
[![Dev](https://img.shields.io/badge/docs-dev-blue.svg)](https://Astroshaper.github.io/AsteroidThermoPhysicalModels.jl/dev)
[![Build Status](https://github.com/Astroshaper/AsteroidThermoPhysicalModels.jl/workflows/CI/badge.svg)](https://github.com/Astroshaper/AsteroidThermoPhysicalModels.jl/actions?query=workflow%3ACI+branch%3Amain)
[![codecov](https://codecov.io/gh/Astroshaper/AsteroidThermoPhysicalModels.jl/branch/main/graph/badge.svg?token=dJBiR91dCD)](https://codecov.io/gh/Astroshaper/AsteroidThermoPhysicalModels.jl)
[![Aqua QA](https://raw.githubusercontent.com/JuliaTesting/Aqua.jl/master/badge.svg)](https://github.com/JuliaTesting/Aqua.jl)

`AsteroidThermoPhysicalModels.jl` is a comprehensive Julia-based toolkit for thermophysical modeling (TPM) of asteroids. It allows you to simulate the temperature distribution of asteroids and predict non-gravitational perturbations on their dynamics (Yarkovsky and YORP effects).

## 📚 Documentation

For detailed documentation, please visit:
- [Stable Documentation](https://Astroshaper.github.io/AsteroidThermoPhysicalModels.jl/stable)
- [Development Documentation](https://Astroshaper.github.io/AsteroidThermoPhysicalModels.jl/dev)

Sample notebooks are available in [Astroshaper-examples](https://github.com/Astroshaper/Astroshaper-examples).

## 🚀 Installation

```julia
using Pkg
Pkg.add("AsteroidThermoPhysicalModels")
using AsteroidThermoPhysicalModels
```

Or in the Julia REPL package mode:
```
julia> ]  # Press ] to enter package mode
pkg> add AsteroidThermoPhysicalModels
```

## 🔍 Features

### Thermophysical Processes
- **Heat Conduction**: 1-dimensional heat diffusion in depth direction
  - Multiple numerical solvers available (explicit Euler, implicit Euler, and Crank-Nicolson methods)
- **Self-Shadowing**: Local shadows cast by topography
- **Self-Heating**: Re-absorption of scattered and radiated photons by surrounding facets
- **Binary Systems**: Support for mutual shadowing (eclipses) and mutual heating between primary and secondary bodies

### Shape Models
- Supports Wavefront OBJ format (*.obj)

### Non-Gravitational Effects
- **Yarkovsky Effect**: Orbital perturbation due to asymmetric thermal emission
- **YORP Effect**: Rotational perturbation due to asymmetric thermal emission

## 🆕 What's New in v0.2.x

### v0.2.1

- `*Ephemerides` constructors now accept `AbstractRange` for `times` — no need to call `collect` beforehand
- `*Ephemerides` constructors now auto-convert plain arrays to `SVector`/`SMatrix` — no need to import `StaticArrays`

### v0.2.0

This release introduces a **Problem-Solver API** inspired by `DifferentialEquations.jl`, replacing the old `run_TPM!` function.

### API Redesign (Breaking)

```julia
# v0.1.x
ephem  = (time = times, sun = r_sun)   # NamedTuple
stpm   = SingleAsteroidThermoPhysicalModel(shape, thermo_params; ...)
result = run_TPM!(stpm, ephem, times_to_save, face_ID)

# v0.2.0
ephem   = SingleAsteroidEphemerides(times, r_sun)
problem = SingleAsteroidThermoPhysicalProblem(shape, thermo_params; ...)
output  = SingleAsteroidOutputSpec(output_times, subsurface_face_ids;
    save_surface_temperature    = true,
    save_subsurface_temperature = true,
    save_face_forces            = false,
    save_forces                 = false,  # true requires R_body_to_inertial in ephem
    save_torques                = false,
)
solution = solve(problem, CrankNicolson();
    ephem               = ephem,
    output              = output,
    initial_temperature = 200.0,
)
```

See the [Migration Guide](CHANGELOG.md#migration-guide) and [CHANGELOG](CHANGELOG.md) for details.

## 🌟 Example

Temperature distribution of asteroid Didymos and its satellite Dimorphos:

<p align="center">
  <img src="https://github.com/user-attachments/assets/a8fd7ad3-5722-4b9b-839d-5c8c2eba6cad" alt="TPM_Didymos">
</p>

## 📖 Basic Usage

The workflow follows a **Problem → Solve → Export** pattern.

```julia
using AsteroidShapeModels
using AsteroidThermoPhysicalModels

# --- Shape model ---
shape = load_shape_obj("path/to/shape.obj"; scale=1000, with_face_visibility=true, with_bvh=true)

# --- Ephemerides ---
# `times`  : epochs [s], any AbstractRange or Vector{Float64}
# `r_sun`  : Sun position in body-fixed frame [m], Vector of length-3 arrays
#   (typically computed from SPICE kernels — see integration test examples)
ephem = SingleAsteroidEphemerides(times, r_sun)

# --- Thermal parameters ---
k     = 0.1    # Thermal conductivity [W/m/K]
ρ     = 1270.0 # Density [kg/m³]
Cₚ    = 600.0  # Heat capacity [J/kg/K]
R_vis = 0.04   # Reflectance in visible light [-]
R_ir  = 0.0    # Reflectance in thermal infrared [-]
ε     = 1.0    # Emissivity [-]
z_max   = 0.6  # Lower boundary depth [m]
n_depth = 61   # Number of depth nodes
Δz      = z_max / (n_depth - 1)

thermo_params = ThermoParams(k, ρ, Cₚ, R_vis, R_ir, ε, z_max, Δz, n_depth)

# --- Problem definition ---
problem = SingleAsteroidThermoPhysicalProblem(shape, thermo_params;
    with_self_shadowing      = true,
    with_self_heating        = true,
    upper_boundary_condition = RadiationBoundaryCondition(),
    lower_boundary_condition = InsulationBoundaryCondition(),
)

# --- Output specification ---
output_times        = ephem.times[end-119:end]  # final rotation period
subsurface_face_ids = [1, 2, 3]                 # faces for saving subsurface temperature profiles
output = SingleAsteroidOutputSpec(output_times, subsurface_face_ids)

# --- Solve ---
solution = solve(problem, CrankNicolson();
    ephem               = ephem,
    output              = output,
    initial_temperature = 200.0,  # [K]
)

# --- Export results ---
export_solution("output/", solution)
# Writes: diagnostics.csv, surface_temperature.csv, subsurface_temperature.csv
```

For a complete end-to-end example with SPICE ephemerides, see [test/TPM_Ryugu/TPM_Ryugu.jl](test/TPM_Ryugu/TPM_Ryugu.jl).

## 📊 Output

The package produces detailed output files including:
- Surface and subsurface temperature distributions
- Thermal forces and torques
- Energy conservation metrics

## 🤝 Contributing

Contributions are welcome! Please feel free to submit a Pull Request.
