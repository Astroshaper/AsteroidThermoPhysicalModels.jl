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
N-body integrators will be implemeted.

## Spin dynamics
Thermal recoil torque (i.e. YORP effect) can be calculated based on a 3-dimensional shape of an asteroid.

## Gravity calculation for asteoids
Gravity filed of a constant-density polyhedron.
