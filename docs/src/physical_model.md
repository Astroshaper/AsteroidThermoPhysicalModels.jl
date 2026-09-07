# Physical Model

## Overview

`AsteroidThermoPhysicalModels.jl` is a comprehensive toolkit for thermophysical modeling of asteroids. This package allows you to simulate the temperature distribution on an asteroid, and predict the non-gravitational force (i.e., Yarkovsky and YORP effects).

The thermophysical model (TPM) considers the following physical processes:

1. **Heat Conduction**: Solves a one-dimensional heat conduction equation to model heat transfer from the surface into the interior of the asteroid.
2. **Self-Shadowing**: Accounts for local shadows cast by topography.
3. **Self-Heating**: Considers re-absorption of scattered light and thermal radiation from surrounding surfaces.
4. **Mutual Shadowing**: For a binary asteroid, accounts for eclipses between the primary and secondary bodies.
5. **Mutual Heating**: For a binary asteroid, considers thermal exchange between the primary and secondary bodies.

## Symbols

The following symbols are used throughout the package to represent various physical quantities:

| Symbol | Unit | Description |
| :----- | :--- | :---------- |
| ``t``                | ``[\mathrm{s}]``                 | Time |
| ``T``                | ``[\mathrm{K}]``                 | Temperature |
| ``R_\text{vis}``     | ``[\text{-}]``                   | Reflectance for visible light |
| ``R_\text{ir}``      | ``[\text{-}]``                   | Reflectance for thermal infrared |
| ``F_\text{sun}``     | ``[\mathrm{W/m^2}]``             | Flux of direct sunlight |
| ``F_\text{scat}``    | ``[\mathrm{W/m^2}]``             | Flux of scattered light |
| ``F_\text{rad}``     | ``[\mathrm{W/m^2}]``             | Flux of thermal radiation from surrounding surface |
| ``\rho``             | ``[\mathrm{kg/m^3}]``            | Density |
| ``C_p``              | ``[\mathrm{J/K}]``               | Heat capacity at constant pressure |
| ``P``                | ``[\mathrm{s}]``                 | Rotation period |
| ``l``                | ``[\mathrm{m}]``                 | Thermal skin depth |
| ``k``                | ``[\mathrm{W/(m \cdot K)}]``     | Thermal conductivity |
| ``z``                | ``[\mathrm{m}]``                 | Depth |
| ``E``                | ``[\mathrm{J}]``                 | Emittance energy |
| ``\Gamma``           | ``[\mathrm{tiu}] = [\mathrm{J \cdot m^{-2} \cdot K^{-1} \cdot s^{-1/2}}]`` | Thermal inertia (cf. [Thermal inertia SI unit proposal](https://nathaniel.putzig.com/research/tiu.html))    |
| ``\varepsilon``      | ``[\text{-}]``                   | Emissivity |
| ``\Phi``             | ``[\mathrm{W/m^2}]``             | Solar energy flux |
| ``\sigma_\text{SB}`` | ``[\mathrm{W/(m^2 \cdot K^4)}]`` | Stefan-Boltzmann constant |

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

### Numerical Solvers

`AsteroidThermoPhysicalModels.jl` provides three numerical methods to solve the heat conduction equation:

1. **Explicit Euler Method** (`ExplicitEulerSolver`)
   - Forward difference in time
   - Conditionally stable: requires ``\lambda = \alpha \Delta t / \Delta z^2 < 0.5``
   - First-order accurate in time
   - Fast for small time steps

2. **Implicit Euler Method** (`ImplicitEulerSolver`)
   - Backward difference in time
   - Unconditionally stable for any time step
   - First-order accurate in time
   - Requires solving a tridiagonal system

3. **Crank-Nicolson Method** (`CrankNicolsonSolver`)
   - Average of forward and backward differences
   - Unconditionally stable
   - Second-order accurate in both time and space
   - Best balance of accuracy and stability

The solver can be specified when creating the thermophysical model:

```julia
# Example: Using Crank-Nicolson solver
stpm = SingleAsteroidTPM(shape, thermo_params;
   SELF_SHADOWING = true,
   SELF_HEATING   = true,
   SOLVER         = AsteroidThermoPhysicalModels.CrankNicolsonSolver(thermo_params),
   BC_UPPER       = AsteroidThermoPhysicalModels.RadiationBoundaryCondition(),
   BC_LOWER       = AsteroidThermoPhysicalModels.InsulationBoundaryCondition(),
)
```

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

### Net force and torque

The photon recoil on each facet ``i`` is integrated over the surface as

```math
\mathbf{F} = \sum_i \mathbf{F}_i, \qquad
\boldsymbol{\tau} = \sum_i \mathbf{r}_i \times \mathbf{F}_i
```

where ``\mathbf{F}_i`` is the recoil force on facet ``i`` (direct emission and reflection, plus the momentum of photons intercepted by other facets when self-heating is enabled) and ``\mathbf{r}_i`` is the position of the facet centre. The net force is the plain vector sum: where a force acts does not enter the motion of the centre of mass, only the torque. The torque is taken about the body-fixed origin, which is therefore assumed to coincide with the centre of mass.

!!! warning "Net force before v0.3.0"
    Up to v0.2.1 the net force was accumulated as ``\sum_i (\hat{\mathbf{r}}_i \cdot \mathbf{F}_i)\,\hat{\mathbf{r}}_i`` — each facet force projected onto the direction of its centre from the origin. This has no physical basis and dropped the tangential part of every facet force. It is harmless on a sphere centred at the origin, where ``\hat{\mathbf{r}}_i = \hat{\mathbf{n}}_i``, and on symmetric polyhedra whose facet centres lie along their normals, but on an irregular body it biases the Yarkovsky force: on the 49k-facet Ryugu shape the rotation-averaged net force came out about 7 % too small and 6° off in direction, and its component along the spin axis had the wrong sign. Net forces computed with earlier versions on non-spherical shapes should be recomputed.

## Binary Asteroid Systems

For binary asteroid systems, `AsteroidThermoPhysicalModels.jl` provides comprehensive modeling of thermal interactions between the primary and secondary bodies.

### Coordinate Systems

The package uses the following coordinate conventions:
- **Primary-fixed frame**: The reference frame fixed to the primary body
- **Secondary-fixed frame**: The reference frame fixed to the secondary body
- **r₁₂**: Position vector from primary to secondary in the primary-fixed frame
- **R₁₂**: Rotation matrix from primary to secondary frame
- **R₂₁ = R₁₂ᵀ**: Rotation matrix from secondary to primary frame

### Coordinate Transformations

For binary systems, coordinate transformations are handled automatically by the unified API:
```julia
# The unified API handles all transformations internally
update_flux_all!(btpm, r☉₁, r₁₂, R₁₂)
```

The package internally computes:
- Sun position in secondary frame: `r☉₂ = R₁₂ * (r☉₁ - r₁₂)`
- Primary position in secondary frame: `r₂₁ = -R₁₂ * r₁₂`

These transformations ensure accurate calculation of mutual shadowing and heating effects.
