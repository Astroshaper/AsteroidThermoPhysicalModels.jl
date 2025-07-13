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

### Coming Soon (v0.2.0)
- Surface roughness modeling integration with `AsteroidShapeModels.jl`
- Enhanced heat conduction solver validation and benchmarks

## 🌟 Example

Temperature distribution of asteroid Didymos and its satellite Dimorphos:

<p align="center">
  <img src="https://github.com/user-attachments/assets/a8fd7ad3-5722-4b9b-839d-5c8c2eba6cad" alt="TPM_Didymos">
</p>

## 📖 Basic Usage

```julia
using AsteroidShapeModels
using AsteroidThermoPhysicalModels


# Load an asteroid shape model
# - `path/to/shape.obj` is the path to your OBJ file (mandatory)
# - `scale` : scale factor for the shape model (e.g., 1000 for km to m conversion)
shape = load_shape_obj("path/to/shape.obj"; scale=1000)

# Set thermal parameters
thermo_params = ThermoParams(
    8.0 * 3600,  # Rotation period [s]
    0.05,        # Thermal skin depth [m]
    200.0,       # Thermal inertia [J m⁻² K⁻¹ s⁻¹/²]
    0.1,         # Reflectance in visible light [-]
    0.0,         # Reflectance in thermal infrared [-]
    0.9,         # Emissivity [-]
    0.5,         # Depth of lower boundary [m]
    0.0125,      # Depth step width [m]
    41           # Number of depth steps
)

# Create TPM model with solver selection
stpm = SingleAsteroidTPM(shape, thermo_params;
    SELF_SHADOWING = true,
    SELF_HEATING   = true,
    SOLVER         = CrankNicolsonSolver(thermo_params),  # Choose solver
    BC_UPPER       = RadiationBoundaryCondition(),
    BC_LOWER       = InsulationBoundaryCondition()
)

# Available solvers:
# - ExplicitEulerSolver(thermo_params)   # Fast but requires small time steps
# - ImplicitEulerSolver(thermo_params)   # Stable for any time step
# - CrankNicolsonSolver(thermo_params)   # Best accuracy

# Run simulation - see documentation for complete examples
```

## 📊 Output

The package produces detailed output files including:
- Surface and subsurface temperature distributions
- Thermal forces and torques
- Energy conservation metrics

## 🤝 Contributing

Contributions are welcome! Please feel free to submit a Pull Request.
