# AsteroidThermoPhysicalModels.jl

`AsteroidThermoPhysicalModels.jl` is a comprehensive Julia-based toolkit for thermophysical modeling (TPM) of asteroids. It allows you to simulate the temperature distribution of asteroids and predict non-gravitational perturbations on their dynamics. Sample notebooks are available in [Astroshaper-examples](https://github.com/Astroshaper/Astroshaper-examples).

## Getting Started

### Installation

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

### Update and Test

To update the package and run tests:

```julia
Pkg.update("AsteroidThermoPhysicalModels")
Pkg.test("AsteroidThermoPhysicalModels")
```

## TPM Workflow

A thermophysical simulation follows this general workflow:

1. **Set up the environment**:

   - Define ephemerides of the target asteroid(s)
   - Load a shape model(s)
   - Set thermophysical parameters

2. **Configure the TPM model**:

   - Choose a solver for the heat conduction equation
   - Set boundary conditions
   - Initialize temperature distribution

3. **Run the simulation**:

   - Execute the TPM calculation
   - Export results

4. **Analyze the results**:

   - Process temperature distributions
   - Calculate non-gravitational effects

## Basic Usage Example

Here's a simple example of setting up and running a thermophysical model for a single asteroid:

```julia
using AsteroidThermoPhysicalModels

# Load shape model
shape = load_shape_obj("asteroid_shape.obj"; scale=1000, with_face_visibility=true, with_bvh=true)

# Ephemerides: one rotation with a synthetic sun trajectory (use SPICE kernels in practice)
P       = 8.0 * 3600  # Rotation period [s]
n_steps = 73
times   = range(0.0, P; length=n_steps)
r_sun   = [[cos(2π*t/P), sin(2π*t/P), 0.0] * 1.5e11 for t in times]  # Sun position in the body-fixed frame [m]
ephem   = SingleAsteroidEphemerides(collect(times), r_sun)

# Material properties and depth grid
thermo_params = ThermoParams(
    conductivity    = 0.1,     # [W/m/K]
    density         = 1270.0,  # [kg/m³]
    heat_capacity   = 600.0,   # [J/kg/K]
    reflectance_vis = 0.1,     # [-]
    reflectance_ir  = 0.0,     # [-]
    emissivity      = 0.9,     # [-]
)
grid_params = GridParams(; z_max=0.5, n_depth=41)  # lower boundary at 0.5 m, 41 depth nodes

# Problem definition
problem = SingleAsteroidThermoPhysicalProblem(shape, thermo_params, grid_params;
    with_self_shadowing      = true,
    with_self_heating        = true,
    upper_boundary_condition = RadiationBoundaryCondition(),
    lower_boundary_condition = InsulationBoundaryCondition(),
)

# Output specification: save the last 10 epochs; subsurface profiles of faces 1–3
output = SingleAsteroidOutputSpec(times[end-9:end]; subsurface_face_ids=[1, 2, 3])

# Solve
solution = solve(problem, CrankNicolson();
    ephem               = ephem,
    output              = output,
    initial_temperature = 200.0,  # [K]
)

# Export results
export_solution("output/", solution)
```

## Output Files

`export_solution` writes the saved quantities as CSV files. Which files appear depends on
the output specification (`SingleAsteroidOutputSpec` / `BinaryAsteroidOutputSpec`):

- `diagnostics.csv` : Energy balance of every saved epoch (`E_in`, `E_out`, `E_cons`)
- `surface_temperature.csv` : Surface temperature of every face [K] as a function of time
- `subsurface_temperature.csv` : Temperature [K] as a function of depth [m] and time (faces in `subsurface_face_ids`)
- `roughness_surface_temperature.csv` : Sub-facet temperatures of the roughness models (with `save_roughness_temperature = true`)
- `thermal_net_forces.csv` : Net thermal force [N] and torque [N·m] in the inertial frame (with `save_forces` / `save_torques`; requires the rotation matrices to the inertial frame in the ephemerides)
- `thermal_face_forces.csv` : Thermal force on every face [N] in the body-fixed frame (with `save_face_forces = true`)

## Binary Asteroid Systems

`AsteroidThermoPhysicalModels.jl` also supports thermophysical modeling of binary asteroid systems. The package can account for mutual shadowing (eclipses) and mutual heating between the primary and secondary bodies.

Example setup for a binary asteroid system:

```julia
# Problems for the primary and secondary
prob1 = SingleAsteroidThermoPhysicalProblem(shape1, thermo_params1, grid_params1; ...)
prob2 = SingleAsteroidThermoPhysicalProblem(shape2, thermo_params2, grid_params2; ...)

# Binary problem
problem = BinaryAsteroidThermoPhysicalProblem(prob1, prob2;
    with_mutual_shadowing = true,
    with_mutual_heating   = true,
)

# Solve (binary_ephem::BinaryAsteroidEphemerides carries the mutual orbit and rotations)
solution = solve(problem, CrankNicolson();
    ephem                         = binary_ephem,
    output                        = BinaryAsteroidOutputSpec(output_times, output_times),
    initial_temperature_primary   = 200.0,
    initial_temperature_secondary = 200.0,
)
```

## Non-Gravitational Effects

### Yarkovsky Effect

The Yarkovsky effect is an orbital perturbation caused by the asymmetric thermal radiation resulting from thermal inertia. This effect primarily affects the semi-major axis of the asteroid's orbit. The thermal force calculated by the TPM can be used to determine the magnitude and direction of the Yarkovsky effect on the asteroid's orbit. You should note that the output force is calculated in the asteroid-fixed frame, so you need to transform it to the inertial frame.

### YORP Effect

The YORP effect (Yarkovsky–O'Keefe–Radzievskii–Paddack effect) is a rotational perturbation resulting from thermal radiation due to the asymmetric shape of the asteroid. This effect influences the rotation rate and the orientation of the asteroid's spin axis. The thermal torque calculated by the TPM can be used to determine the magnitude and direction of the YORP effect on the asteroid's rotation. You should note that the output torque is calculated in the asteroid-fixed frame, so you need to transform it to the inertial frame.

## References

For more details on the package development and the thermophysical model, refer to the following references:

- [Kanamaru et al. (2024)](https://doi.org/10.57350/jesa.206) - Thermophysical Model Development for Hera Mission to Simulate Non-Gravitational Acceleration on Binary Asteroid
- [金丸ほか (2024)](https://www.wakusei.jp/book/pp/2024/2024-3/235-242_p.pdf) - 二重小惑星探査計画Heraに向けた小惑星熱物理モデルの開発 (in Japanese)

The package contributes the following publications:
- [Kanamaru et al. (2021)](https://doi.org/10.1029/2021JE006863) - YORP Effect on Asteroid 162173 Ryugu: Implications for the Dynamical History
- [Zhou et al. (2024)](https://iopscience.iop.org/article/10.3847/2041-8213/ad4f7f) - The Yarkovsky Effect on the Long-term Evolution of Binary Asteroids
