# Physical Model

## Overview

`AsteroidThermoPhysicalModels.jl` is a comprehensive toolkit for thermophysical modeling of asteroids. This package allows you to simulate the temperature distribution on an asteroid, and predict the non-gravitational force (i.e., Yarkovsky and YORP effects).

The thermophysical model (TPM) considers the following physical processes:

1. **Heat Conduction**: Solves a one-dimensional heat conduction equation to model heat transfer from the surface into the interior of the asteroid.
2. **Self-Shadowing**: Accounts for local shadows cast by topography.
3. **Self-Heating**: Considers re-absorption of scattered light and thermal radiation from surrounding surfaces.
4. **Mutual Shadowing**: For a binary asteroid, accounts for eclipses between the primary and secondary bodies.
5. **Mutual Heating**: For a binary asteroid, considers thermal exchange between the primary and secondary bodies.

## Heat Conduction Equation

Heat conduction within the asteroid is modeled by the following one-dimensional heat diffusion equation:

```math
\rho C_p \frac{\partial T}{\partial t} = k \frac{\partial^2 T}{\partial z^2}
```

where:
- ``T(z,t)`` is the temperature at depth ``z`` and time ``t``
- ``\rho`` is the density
- ``C_p`` is the specific heat capacity at constant pressure
- ``k`` is the thermal conductivity

### Boundary Conditions

#### Upper Boundary Condition (Surface)

At the surface, a radiative equilibrium boundary condition is applied:

```math
-k \frac{\partial T}{\partial z}\bigg|_{z=0} = (1-R_\text{vis})(F_\text{sun} + F_\text{scat}) + (1-R_\text{ir})F_\text{rad} - \varepsilon \sigma T^4
```

where:
- The left side represents the heat flux from the surface into the interior
- The right side represents the energy balance at the surface (absorbed solar radiation, scattered light, and thermal radiation minus emitted thermal radiation)

#### Lower Boundary Condition

At the lower boundary, an insulation condition is typically applied:

```math
\frac{\partial T}{\partial z}\bigg|_{z=z_\text{max}} = 0
```

## Thermal Inertia

Thermal inertia is a physical quantity that represents the ability of a material to resist temperature changes, defined by:

```math
\Gamma = \sqrt{k \rho C_p}
```

The unit is `tiu` (thermal inertia unit) or `J·m⁻²·K⁻¹·s⁻¹/²`.

## Thermal Skin Depth

The thermal skin depth is a characteristic length that represents how far a periodic thermal wave penetrates into a material:

```math
l = \sqrt{\frac{4\pi P k}{\rho C_p}}
```

where ``P`` is the period of the thermal cycle (typically the rotation period of the asteroid).

## Non-Gravitational Effects

### Yarkovsky Effect

The Yarkovsky effect is an orbital perturbation caused by the asymmetric thermal emission resulting from the day-night temperature difference due to the asteroid's rotation. This effect primarily affects the semi-major axis of the asteroid's orbit.

### YORP Effect

The YORP effect (Yarkovsky-O'Keefe-Radzievskii-Paddack effect) is a rotational perturbation resulting from thermal emission due to the asymmetric shape of the asteroid. This effect influences the rotation rate and the orientation of the asteroid's spin axis.

## Symbols

| Symbol | Unit | Description |
| :----- | :--- | :---------- |
| ``t`` | ``[\mathrm{s}]`` | Time |
| ``T`` | ``[\mathrm{K}]`` | Temperature |
| ``R_\text{vis}`` | ``[\text{-}]`` | Reflectance for visible light |
| ``R_\text{ir}`` | ``[\text{-}]`` | Reflectance for thermal infrared |
| ``F_\text{sun}`` | ``[\mathrm{W/m^2}]`` | Flux of direct sunlight |
| ``F_\text{scat}`` | ``[\mathrm{W/m^2}]`` | Flux of scattered light |
| ``F_\text{rad}`` | ``[\mathrm{W/m^2}]`` | Flux of thermal radiation from surrounding surface |
| ``\rho`` | ``[\mathrm{kg/m^3}]`` | Density |
| ``C_p`` | ``[\mathrm{J/(kg \cdot K)}]`` | Heat capacity at constant pressure |
| ``P`` | ``[\mathrm{s}]`` | Rotation period |
| ``l`` | ``[\mathrm{m}]`` | Thermal skin depth |
| ``k`` | ``[\mathrm{W/(m \cdot K)}]`` | Thermal conductivity |
| ``z`` | ``[\mathrm{m}]`` | Depth |
| ``E`` | ``[\mathrm{J}]`` | Emittance energy |
| ``\Gamma`` | ``[\mathrm{tiu}] = [\mathrm{J \cdot m^{-2} \cdot K^{-1} \cdot s^{-1/2}}]`` | Thermal inertia (cf. [Thermal inertia SI unit proposal](https://nathaniel.putzig.com/research/tiu.html)) |
| ``\varepsilon`` | ``[\text{-}]`` | Emissivity |
| ``\Phi`` | ``[\mathrm{W/m^2}]`` | Solar energy flux |
| ``\sigma_\text{SB}`` | ``[\mathrm{W/(m^2 \cdot K^4)}]`` | Stefan-Boltzmann constant |
