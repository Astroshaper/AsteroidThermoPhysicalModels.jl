
![test](https://github.com/MasanoriKanamaru/Astroshaper/files/7840328/Astroshaper_logo.pdf)


# Astroshaper

[![Build Status](https://travis-ci.com/MasanoriKanamaru/Astroshaper.svg?branch=main)](https://travis-ci.com/MasanoriKanamaru/Astroshaper)

Julia-based toolkit for dynamical simulations of planets and small solar system bodies. Tutorials that you can easily try will be prepared soon.

## Installation

    using Pkg
    Pkg.add(url="https://github.com/MasanoriKanamaru/Astroshaper")
    using Astroshaper

You can update the module and run tests as follows.

    Pkg.update("Astroshaper")
    Pkg.test("Astroshaper")

## Orbital dynamics
You can simulate orbital evolution of planets and small bodies under gravity interaction and various perturbations.
As for the orbital integrators, you can choose from Euler, leapfrog,  4th-degree Hermite methods (Note that my implementation of the Hermite method is being verified). Thermophysical perturbation on orbital motion of an asteroid, such as Yarkovsky effect (Bottke et al., 2006), will be implemented.


## Spin dynamics
Based on the thermophysics of an airless rocky body, you can simulate the distribution of the surface temperature and thermal recoil torque on the body, i.e. YORP effect (Rubincam, 2000; Bottke et al., 2006). The simulation can receive a 3-dimensional shape model of an asteroid (or local elevation model) in the Wavefront OBJ format (*.obj).

The code currently includes the thermophysical effects as follows:
- 1-dimensional heat diffusion in depth direction
- Self-shadowing: Detection of facets blocking solar rays
- Self-heating: Re-absorption of scatterd and radiated photons by surrounding facets. Only sigle scattering is implemented.

## Gravity calculation for asteoids
You can calculate the precise gravity field of an irregularly shaped body, based on the constant-density polyhedron method (Werner & Scheeres, 1997).

## Start to play
Let's visualize a shape model of asteroid Ryugu.
Please downlad a Ryugu model from Astroshaper/test/ryugu_test.obj.

    using Astroshaper

    shapepath = "ryugu_test.obj"                  # Path to the shape model
    shape = Shape(shapepath; scale=1000, find_visible_facets=true)

    draw(shape)
    draw(shape, data=:radius)                     # Radius of each surface facet
    draw(shape; data=:illumination, r̂☉=[1,0,0.])  # Illumination when the Sun is in the direction of r̂☉

<img width="853" alt="start_to_play" src="https://user-images.githubusercontent.com/21192162/148867718-3b3bc681-c210-4631-9c64-263be16161a7.png">
