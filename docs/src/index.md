# AsteroidThermoPhysicalModels.jl

`AsteroidThermoPhysicalModels.jl` is a comprehensive Julia-based toolkit for thermophysical modeling (TPM) of asteroids. It allows you to simulate the temperature distribution of asteroids and predict non-gravitational perturbations on their dynamics. Sample notebooks are available in [Astroshaper-examples](https://github.com/Astroshaper/Astroshaper-examples).

## Getting Started

### Installation

To install `AsteroidThermoPhysicalModels.jl`, use Julia's package manager:

```julia
using Pkg
Pkg.add(url="https://github.com/Astroshaper/AsteroidThermoPhysicalModels.jl")
```

If you'd like to install a specific version, you can add the version number:

```julia
using Pkg
Pkg.add(url="https://github.com/Astroshaper/AsteroidThermoPhysicalModels.jl#0.0.6") 
```

Then, load the package:
```julia
using AsteroidThermoPhysicalModels
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
shape = load_shape_obj("asteroid_shape.obj"; scale=1000, find_visible_facets=true)

# Set thermal parameters
P = 8.0 * 3600  # Rotation period [s]
l = 0.05         # Thermal skin depth [m]
Γ = 200.0        # Thermal inertia [J m⁻² K⁻¹ s⁻¹/²]
R_vis = 0.1      # Reflectance in visible light [-]
R_ir = 0.0       # Reflectance in thermal infrared [-]
ε = 0.9          # Emissivity [-]
z_max = 0.5      # Depth of lower boundary [m]
n_depth = 41     # Number of depth steps
Δz = z_max / (n_depth - 1)  # Depth step width [m]

thermo_params = ThermoParams(P, l, Γ, R_vis, R_ir, ε, z_max, Δz, n_depth)

# Create TPM model
stpm = SingleAsteroidThermoPhysicalModel(shape, thermo_params;
    SELF_SHADOWING = true,
    SELF_HEATING = true,
    SOLVER = ForwardEulerSolver(thermo_params),
    BC_UPPER = RadiationBoundaryCondition(),
    BC_LOWER = InsulationBoundaryCondition()
)

# Initialize temperature
init_temperature!(stpm, 200.0)  # Initial temperature [K]

# Define ephemerides (simplified example)
n_steps = 72
times = range(0.0, P; length=n_steps)
sun_positions = [normalize([cos(2π*t/P), sin(2π*t/P), 0.0]) * 1.5e11 for t in times]

ephem = (time = collect(times), sun = sun_positions)

# Define which timesteps to save detailed data for
times_to_save = times[end-10:end]  # Save last 10 timesteps
face_ID = [1, 2, 3]  # Save subsurface temperature for these faces

# Run TPM
result = run_TPM!(stpm, ephem, times_to_save, face_ID)

# Export results
export_TPM_results("output_directory", result)
```

## Output Files

After running the simulation, you will get the following output files. Note that the thermal force and torque are calculated in the asteroid-fixed frame.

- `physical_quantities.csv`
    - `time`     : Time steps [s]
    - `E_in`     : Input energy per second on the whole surface [W]
    - `E_out`    : Output energy per second from the whole surface [W]
    - `E_cons`   : Energy conservation ratio [-], ratio of total energy going out to total energy coming in during the last rotation cycle
    - `force_x`  : x-component of the thermal force [N]
    - `force_y`  : y-component of the thermal force [N]
    - `force_z`  : z-component of the thermal force [N]
    - `torque_x` : x-component of the thermal torque [N·m]
    - `torque_y` : y-component of the thermal torque [N·m]
    - `torque_z` : z-component of the thermal torque [N·m]

- `subsurface_temperature.csv` : Temperature [K] as a function of depth [m] and time [s]
- `surface_temperature.csv` : Surface temperature of every face [K] as a function of time [s]
- `thermal_force.csv` : Thermal force on every face of the shape model [N] as a function of time

## Binary Asteroid Systems

`AsteroidThermoPhysicalModels.jl` also supports thermophysical modeling of binary asteroid systems. The package can account for mutual shadowing (eclipses) and mutual heating between the primary and secondary bodies.

Example setup for a binary asteroid system:

```julia
# Create TPM models for primary and secondary
primary_tpm = SingleAsteroidThermoPhysicalModel(shape_primary, thermo_params_primary; ...)
secondary_tpm = SingleAsteroidThermoPhysicalModel(shape_secondary, thermo_params_secondary; ...)

# Create binary TPM model
btpm = BinaryAsteroidThermoPhysicalModel(primary_tpm, secondary_tpm;
    MUTUAL_SHADOWING = true,
    MUTUAL_HEATING = true
)

# Run TPM for binary system
result = run_TPM!(btpm, binary_ephem, times_to_save, face_ID_pri, face_ID_sec)
```

## Non-Gravitational Effects

### Yarkovsky Effect

The Yarkovsky effect is an orbital perturbation caused by the asymmetric thermal radiation resulting from thermal inertia. This effect primarily affects the semi-major axis of the asteroid's orbit. The thermal force calculated by the TPM can be used to determine the magnitude and direction of the Yarkovsky effect on the asteroid's orbit. You should note that the output force is calculated in the asteroid-fixed frame, so you need to transform it to the inertial frame.

### YORP Effect

The YORP effect (Yarkovsky–O'Keefe–Radzievskii–Paddack effect) is a rotational perturbation resulting from thermal radiatoion due to the asymmetric shape of the asteroid. This effect influences the rotation rate and the orientation of the asteroid's spin axis. The thermal torque calculated by the TPM can be used to determine the magnitude and direction of the YORP effect on the asteroid's rotation. You should note that the output torque is calculated in the asteroid-fixed frame, so you need to transform it to the inertial frame.

## References

For more details on the package development and the thermophysical model, refer to the following references:

- [Kanamaru et al. (2024)](https://doi.org/10.57350/jesa.206) - Thermophysical Model Development for Hera Mission to Simulate Non-Gravitational Acceleration on Binary Asteroid
- [金丸ほか (2024)](https://www.wakusei.jp/book/pp/2024/2024-3/235-242_p.pdf) - 二重小惑星探査計画Heraに向けた小惑星熱物理モデルの開発 (in Japanese)

The package contributes the following publications:
- [Kanamaru et al. (2021)](https://doi.org/10.1029/2021JE006863) - YORP Effect on Asteroid 162173 Ryugu: Implications for the Dynamical History
- [Zhou et al. (2024)](https://iopscience.iop.org/article/10.3847/2041-8213/ad4f7f) - The Yarkovsky Effect on the Long-term Evolution of Binary Asteroids
