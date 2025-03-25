# AsteroidThermoPhysicalModels.jl

[![Stable](https://img.shields.io/badge/docs-stable-blue.svg)](https://Astroshaper.github.io/AsteroidThermoPhysicalModels.jl/stable)
[![Dev](https://img.shields.io/badge/docs-dev-blue.svg)](https://Astroshaper.github.io/AsteroidThermoPhysicalModels.jl/dev)
[![Build Status](https://github.com/Astroshaper/AsteroidThermoPhysicalModels.jl/workflows/CI/badge.svg)](https://github.com/Astroshaper/AsteroidThermoPhysicalModels.jl/actions?query=workflow%3ACI+branch%3Amain)
[![codecov](https://codecov.io/gh/Astroshaper/AsteroidThermoPhysicalModels.jl/branch/main/graph/badge.svg?token=dJBiR91dCD)](https://codecov.io/gh/Astroshaper/AsteroidThermoPhysicalModels.jl)
[![Aqua QA](https://raw.githubusercontent.com/JuliaTesting/Aqua.jl/master/badge.svg)](https://github.com/JuliaTesting/Aqua.jl)

A Julia-based toolkit for thermophysical modeling of asteroids. It allows you to simulate the temperature distribution of the asteroid and predict non-gravitational perturbations on its dynamics. Sample notebooks are available in [Astroshaper-examples](https://github.com/Astroshaper/Astroshaper-examples).

## Installation

```julia
using Pkg
Pkg.add(url="https://github.com/Astroshaper/AsteroidThermoPhysicalModels.jl")
using AsteroidThermoPhysicalModels
```

You can update the module and run tests as follows.

```julia
Pkg.update("AsteroidThermoPhysicalModels")
Pkg.test("AsteroidThermoPhysicalModels")
```

## Thermophysical modeling
Based on orbit, spin, and 3-D shape, you can calculate the distribution of the surface temperature on an asteroid. The temperature distribution can be used to calculate the non-gravitational perturbations on its orbital and rotational motion (Yarkovsky and YORP effects, respectively).

### Available format for shape model
- Wavefront OBJ format (\*.obj)

### Thermophysics included
- 1-dimensional heat diffusion in depth direction
- Self-shadowing: Local shadows casted by topography
- Self-heating: Re-absorption of scattered and radiated photons by surrounding facets. Only single scattering is implemented.
- Mutual shadowing (eclipse) for a binary asteroid
- Mutual heating for a binary asteroid

### Thermophysics to be implemented
- Surface roughness (smaller than facets of the shape model)

### Example
Temperature distribution of asteroid Didymos and its satellite Dimorphos.

![TPM_Didymos](https://github.com/user-attachments/assets/a8fd7ad3-5722-4b9b-839d-5c8c2eba6cad)

