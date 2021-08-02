# Astroshaper

[![Build Status](https://travis-ci.com/MasanoriKanamaru/Astroshaper.svg?branch=main)](https://travis-ci.com/MasanoriKanamaru/Astroshaper)

Julia-based toolkit for dynamical simulations of planets and small solar system bodies.

## Installation

    using Pkg
    Pkg.add("https://github.com/MasanoriKanamaru/Astroshaper")
    using Astroshaper

You can update the module and run tests as follows.

    Pkg.update("Astroshaper")
    Pkg.test("Astroshaper")

## Orbital dynamics
You can simulate orbital evolution of planets and small bodies under gravity interaction and various perturbations.
As for the orbital integrators, you can choose from Euler, leapfrog,  4th-degree Hermite methods (Note that my implementation of the Hermite method is being verified). Thermophysical perturbation on orbital motion of an asteroid, such as Yarkovsky effect (Bottke et al. 2006), will be implemented.


## Spin dynamics
Thermal recoil torque (i.e. YORP effect) can be calculated based on a 3-dimensional shape of an asteroid.

## Gravity calculation for asteoids
Gravity filed of a constant-density polyhedron.
