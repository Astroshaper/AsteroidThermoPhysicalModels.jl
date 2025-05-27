<div align="center">
  <h1>AsteroidThermoPhysicalModels.jl</h1>
</div>

<p align="center">
  <a href="https://Astroshaper.github.io/AsteroidThermoPhysicalModels.jl/stable">
    <img src="https://img.shields.io/badge/docs-stable-blue.svg" alt="Stable">
  </a>
  <a href="https://Astroshaper.github.io/AsteroidThermoPhysicalModels.jl/dev">
    <img src="https://img.shields.io/badge/docs-dev-blue.svg" alt="Dev">
  </a>
  <a href="https://github.com/Astroshaper/AsteroidThermoPhysicalModels.jl/actions?query=workflow%3ACI+branch%3Amain">
    <img src="https://github.com/Astroshaper/AsteroidThermoPhysicalModels.jl/workflows/CI/badge.svg" alt="Build Status">
  </a>
  <a href="https://codecov.io/gh/Astroshaper/AsteroidThermoPhysicalModels.jl">
    <img src="https://codecov.io/gh/Astroshaper/AsteroidThermoPhysicalModels.jl/branch/main/graph/badge.svg?token=dJBiR91dCD" alt="codecov">
  </a>
  <a href="https://github.com/JuliaTesting/Aqua.jl">
    <img src="https://raw.githubusercontent.com/JuliaTesting/Aqua.jl/master/badge.svg" alt="Aqua QA">
  </a>
</p>

`AsteroidThermoPhysicalModels.jl` is a comprehensive Julia-based toolkit for thermophysical modeling (TPM) of asteroids. It allows you to simulate the temperature distribution of asteroids and predict non-gravitational perturbations on their dynamics (Yarkovsky and YORP effects).

## 📚 Documentation

For detailed documentation, please visit:
- [Stable Documentation](https://Astroshaper.github.io/AsteroidThermoPhysicalModels.jl/stable)
- [Development Documentation](https://Astroshaper.github.io/AsteroidThermoPhysicalModels.jl/dev)

Sample notebooks are available in [Astroshaper-examples](https://github.com/Astroshaper/Astroshaper-examples).

## 🚀 Installation

```julia
using Pkg
Pkg.add(url="https://github.com/Astroshaper/AsteroidThermoPhysicalModels.jl")
using AsteroidThermoPhysicalModels
```

You can update the package and run tests as follows:

```julia
Pkg.update("AsteroidThermoPhysicalModels")
Pkg.test("AsteroidThermoPhysicalModels")
```

## 🔍 Features

### Thermophysical Processes
- **Heat Conduction**: 1-dimensional heat diffusion in depth direction
- **Self-Shadowing**: Local shadows cast by topography
- **Self-Heating**: Re-absorption of scattered and radiated photons by surrounding facets
- **Binary Systems**: Support for mutual shadowing (eclipses) and mutual heating between primary and secondary bodies

### Shape Models
- Supports Wavefront OBJ format (*.obj)

### Non-Gravitational Effects
- **Yarkovsky Effect**: Orbital perturbation due to asymmetric thermal emission
- **YORP Effect**: Rotational perturbation due to asymmetric thermal emission

### Coming Soon
- Surface roughness modeling (at scales smaller than facets of the shape model)

## 🌟 Example

Temperature distribution of asteroid Didymos and its satellite Dimorphos:

<p align="center">
  <img src="https://github.com/user-attachments/assets/a8fd7ad3-5722-4b9b-839d-5c8c2eba6cad" alt="TPM_Didymos">
</p>

## 📖 Basic Usage

```julia
using AsteroidThermoPhysicalModels

# Load shape model
shape = load_shape_obj("asteroid_shape.obj"; scale=1000, find_visible_facets=true)

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

# Create and run TPM model
# See documentation for complete examples
```

## 📊 Output

The package produces detailed output files including:
- Surface and subsurface temperature distributions
- Thermal forces and torques
- Energy conservation metrics

## 🤝 Contributing

Contributions are welcome! Please feel free to submit a Pull Request.
