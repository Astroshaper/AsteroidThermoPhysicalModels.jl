# AsteroidThermoPhysicalModels.jl

A Julia-based toolkit for thermophysical modeling (TPM) of asteroids. It allows you to simulate the temperature distribution of the asteroid and predict non-gravitational perturbations on its dynamics. Sample notebooks are available in [Astroshaper-examples](https://github.com/Astroshaper/Astroshaper-examples).

## TPM workflow
A thermophysical simulation is performed as the following workflow.

- Set ephemerides of a target asteroid(s)
- Load a shape model(s)
- Set thermophysical parameters

After running the simulation, you will get the following output files. You should note that the thermal force and torque is calculated in the asteroid-fixed frame.
- `physical_quantities.csv`
    - `time`     : Time steps
    - `E_in`     : Input energy per second on the whole surface [W]
    - `E_out`    : Output enegey per second from the whole surface [W]
    - `E_cons`   : Energy conservation ratio [-], ratio of total energy going out to total energy coming in during the last rotation cycle
    - `force_x`  : x-component of the thermal force
    - `force_y`  : y-component of the thermal force
    - `force_z`  : z-component of the thermal force
    - `torque_x` : x-component of the thermal torque
    - `torque_y` : y-component of the thermal torque
    - `torque_z` : z-component of the thermal torque
- `subsurface_temperature.csv` : Temperature [K] as a function of depth [m] and time [s]
- `surface_temperature.csv` : Surface temperature of every face [K] as a function of time [s]
- `thermal_force.csv` : Thermal force on every face of the shape model [N] as a function of time

## Surface/sub-surface temperature analysis
Coming soon.

## Non-Gravitational Effects

### Yarkovsky effect
Coming soon.

### YORP effect
Coming soon.

