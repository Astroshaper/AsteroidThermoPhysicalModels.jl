![](logo/Astroshaper_logo.png)

# ThermoPhysicalModeling.jl

[![Stable](https://img.shields.io/badge/docs-stable-blue.svg)](https://MasanoriKanamaru.github.io/ThermoPhysicalModeling.jl/stable)
[![Dev](https://img.shields.io/badge/docs-dev-blue.svg)](https://MasanoriKanamaru.github.io/ThermoPhysicalModeling.jl/dev)
[![Build Status](https://github.com/MasanoriKanamaru/ThermoPhysicalModeling.jl/workflows/CI/badge.svg)](https://github.com/MasanoriKanamaru/ThermoPhysicalModeling.jl/actions?query=workflow%3ACI+branch%3Amain)
[![codecov](https://codecov.io/gh/MasanoriKanamaru/ThermoPhysicalModeling.jl/branch/main/graph/badge.svg?token=dJBiR91dCD)](https://codecov.io/gh/MasanoriKanamaru/ThermoPhysicalModeling.jl)
[![Aqua QA](https://raw.githubusercontent.com/JuliaTesting/Aqua.jl/master/badge.svg)](https://github.com/JuliaTesting/Aqua.jl)

Julia-based toolkit for dynamical simulations of planets and small solar system bodies.

## Installation

    using Pkg
    Pkg.add(url="https://github.com/MasanoriKanamaru/ThermoPhysicalModeling.jl")
    using ThermoPhysicalModeling

You can update the module and run tests as follows.

    Pkg.update("ThermoPhysicalModeling")
    Pkg.test("ThermoPhysicalModeling")

## Orbital dynamics
You can simulate orbital evolution of planets and small bodies under gravity interaction and various perturbations.
As for the orbital integrators, you can choose from Euler, leapfrog,  4th-degree Hermite methods (Note that my implementation of the Hermite method is being verified). Thermophysical perturbation on orbital motion of an asteroid, that is, Yarkovsky effect will be implemented.

### Example
Orbital calculation of the solar system. The gravity of 8 planets and Pluto are considered.

<img src="https://user-images.githubusercontent.com/21192162/149469835-42fed69f-f93e-4123-80ea-baa6497aebca.gif" width="400">

## Themophysical modeling
Based on orbit, spin, and 3-D shape, you can calculate the distribution of the surface temperature on an asteroid. The temperature distribution can be used to calculate the non-gravitational perturbations on its orbital and rotational motion (Yarkovsky and YORP effects, respectively).

### Available format for shape model
- Wavefront OBJ format (\*.obj)

### Thermophysics included
- 1-dimensional heat diffusion in depth direction
- Self-shadowing: Local shadows casted by topography
- Self-heating: Re-absorption of scatterd and radiated photons by surrounding facets. Only sigle scattering is implemented.

### Thermophysics to be implemted
- Surface roughness (smaller than facets of the shape model)

### Example
Distribution of surface temperature on asteroid Ryugu. The color map ranges from 200 to 400 K.

![Thermophysics_Ryugu](https://user-images.githubusercontent.com/21192162/149468024-f403011f-b3d3-47ce-a69c-7daf78a40658.png)


## Gravity calculation for asteoids
You can calculate the precise gravity field of an irregularly shaped body, based on the constant-density polyhedron method (Werner & Scheeres, 1997).

### Example
Distribution of dynamical elevation on asteroid Itokawa. The color map ranges from -25 to 55 m.

![Gravity_Itokawa](https://user-images.githubusercontent.com/21192162/149465150-6cead63e-6027-402f-b866-5111dc5321a7.png)

## Start to play
Let's visualize a shape model of asteroid Ryugu.
Please downlad a Ryugu model from ThermoPhysicalModeling/test/ryugu_test.obj.

    using ThermoPhysicalModeling

    shapepath = "ryugu_test.obj"  # Path to the shape model
    shape = Shape(shapepath; scale=1000, find_visible_facets=true)

    draw(shape)
    # draw(shape, data=:radius)                     # Radius of each surface facet
    # draw(shape; data=:illumination, r̂☉=[1,0,0.])  # Illumination when the Sun is in the direction of r̂☉

<img width="300" alt="start_to_play" src="https://user-images.githubusercontent.com/21192162/148867940-21db4a00-8aef-4030-ab94-397d4f3b572c.png">
