# Usage Examples

This section provides detailed examples of how to use `AsteroidThermoPhysicalModels.jl` for various scenarios.
For more detailed examples, please refer to the [Astroshaper-examples](https://github.com/Astroshaper/Astroshaper-examples) repository.


## Single Asteroid Example (Ryugu)

This example demonstrates how to set up and run a thermophysical model for asteroid Ryugu using SPICE kernels for ephemerides.

```julia
using AsteroidShapeModels
using AsteroidThermoPhysicalModels
using Downloads
using SPICE

##= Download SPICE kernels and shape model =##
paths_kernel = [
    "lsk/naif0012.tls",
    "pck/hyb2_ryugu_shape_v20190328.tpc",
    "fk/hyb2_ryugu_v01.tf",
    "spk/2162173_Ryugu.bsp",
]

for path_kernel in paths_kernel
    url_kernel = "https://data.darts.isas.jaxa.jp/pub/hayabusa2/old/2020/spice_bundle/spice_kernels/$(path_kernel)"
    filepath = joinpath("kernel", path_kernel)
    mkpath(dirname(filepath))
    isfile(filepath) || Downloads.download(url_kernel, filepath)
end

path_shape = "SHAPE_SFM_49k_v20180804.obj"
url_shape  = "https://data.darts.isas.jaxa.jp/pub/hayabusa2/paper/Watanabe_2019/$(path_shape)"
filepath_shape = joinpath("shape", path_shape)
mkpath(dirname(filepath_shape))
isfile(filepath_shape) || Downloads.download(url_shape, filepath_shape)

##= Load SPICE kernels =##
for path_kernel in paths_kernel
    SPICE.furnsh(joinpath("kernel", path_kernel))
end

##= Build ephemerides =##
P = SPICE.convrt(7.63262, "hours", "seconds")  # Rotation period of Ryugu [s]

n_cycle          = 2   # Number of rotation cycles to simulate
n_step_in_cycle  = 72  # Time steps per rotation

et_begin = SPICE.utc2et("2018-07-01T00:00:00")
et_end   = et_begin + P * n_cycle
et_range = range(et_begin, et_end; length=n_step_in_cycle * n_cycle + 1)

r_sun = [SPICE.spkpos("SUN", et, "RYUGU_FIXED", "None", "RYUGU")[1] for et in et_range]
r_sun .*= 1000  # Convert [km] to [m]

ephem = SingleAsteroidEphemerides(et_range, r_sun)

SPICE.kclear()

##= Load shape model =##
shape = load_shape_obj(joinpath("shape", path_shape); scale=1000, with_face_visibility=true, with_bvh=true)

##= Thermal properties and grid settings =##
thermo_params = ThermoParams(
    conductivity    = 0.1,     # Thermal conductivity [W/m/K]
    density         = 1270.0,  # Density [kg/m³]
    heat_capacity   = 600.0,   # Specific heat capacity [J/kg/K]
    reflectance_vis = 0.04,    # Reflectance in visible light [-]
    reflectance_ir  = 0.0,     # Reflectance in thermal infrared [-]
    emissivity      = 1.0,     # Thermal emissivity [-]
)

grid_params = GridParams(;
    z_max   = 0.6,  # Depth of the lower boundary [m]
    n_depth = 61,   # Number of depth nodes
)

##= Define thermophysical problem =##
problem = SingleAsteroidThermoPhysicalProblem(shape, thermo_params, grid_params;
    with_self_shadowing      = true,
    with_self_heating        = true,
    upper_boundary_condition = RadiationBoundaryCondition(),
    lower_boundary_condition = InsulationBoundaryCondition(),
)

##= Output specification =##
output_times        = collect(et_range)[end-n_step_in_cycle:end]  # Save the final rotation
subsurface_face_ids = [1, 2, 3, 4, 10]                            # Face indices for subsurface output

output = SingleAsteroidOutputSpec(output_times, subsurface_face_ids)

##= Run TPM =##
solution = solve(problem, ExplicitEuler();
    ephem               = ephem,
    output              = output,
    initial_temperature = 200.0,
)

export_solution("path/to/save", solution)
```

## Binary Asteroid Example (Didymos-Dimorphos)

This example demonstrates how to set up and run a thermophysical model for the binary asteroid system Didymos-Dimorphos using SPICE kernels for ephemerides.

```julia
using AsteroidShapeModels
using AsteroidThermoPhysicalModels
using Downloads
using SPICE

##= SPICE kernels and shape models =##
paths_kernel = [
    "fk/hera_v10.tf",
    "lsk/naif0012.tls",
    "pck/hera_didymos_v06.tpc",
    "spk/de432s.bsp",
    "spk/didymos_hor_000101_500101_v01.bsp",
    "spk/didymos_gmv_260901_311001_v01.bsp",
]
paths_shape = [
    "g_50677mm_rad_obj_didy_0000n00000_v001.obj",
    "g_08438mm_lgt_obj_dimo_0000n00000_v002.obj",
]

##= Download SPICE kernels =##
for path_kernel in paths_kernel
    url_kernel = "https://s2e2.cosmos.esa.int/bitbucket/projects/SPICE_KERNELS/repos/hera/raw/kernels/$(path_kernel)?at=refs%2Ftags%2Fv161_20230929_001"
    filepath = joinpath("kernel", path_kernel)
    mkpath(dirname(filepath))
    isfile(filepath) || Downloads.download(url_kernel, filepath)
end

##= Download shape models =##
for path_shape in paths_shape
    url_shape = "https://s2e2.cosmos.esa.int/bitbucket/projects/SPICE_KERNELS/repos/hera/raw/kernels/dsk/$(path_shape)?at=refs%2Ftags%2Fv161_20230929_001"
    filepath = joinpath("shape", path_shape)
    mkpath(dirname(filepath))
    isfile(filepath) || Downloads.download(url_shape, filepath)
end

##= Load SPICE kernels =##
for path_kernel in paths_kernel
    SPICE.furnsh(joinpath("kernel", path_kernel))
end

##= Build ephemerides =##
P₂ = SPICE.convrt(11.93, "hours", "seconds")  # Rotation period of Dimorphos [s]

n_cycle         = 2   # Number of rotation cycles to simulate
n_step_in_cycle = 72  # Time steps per rotation

et_begin = SPICE.utc2et("2027-02-18T00:00:00")
et_end   = et_begin + P₂ * n_cycle
et_range = range(et_begin, et_end; length=n_step_in_cycle * n_cycle + 1)

r_sun                  = [SPICE.spkpos("SUN"      , et, "DIDYMOS_FIXED"  , "None", "DIDYMOS")[1] for et in et_range]
r_secondary            = [SPICE.spkpos("DIMORPHOS", et, "DIDYMOS_FIXED"  , "None", "DIDYMOS")[1] for et in et_range]
R_primary_to_secondary = [SPICE.pxform("DIDYMOS_FIXED", "DIMORPHOS_FIXED", et)                          for et in et_range]

r_sun       .*= 1000  # Convert [km] to [m]
r_secondary .*= 1000  # Convert [km] to [m]

ephem = BinaryAsteroidEphemerides(et_range, r_sun, r_secondary, R_primary_to_secondary)

SPICE.kclear()

##= Load shape models =##
shape1 = load_shape_obj(joinpath("shape", paths_shape[1]); scale=1000, with_face_visibility=true, with_bvh=true)
shape2 = load_shape_obj(joinpath("shape", paths_shape[2]); scale=1000, with_face_visibility=true, with_bvh=true)

##= Thermal properties and grid settings =##
# [cf. Michel+2016; Naidu+2020]
thermo_params = ThermoParams(
    conductivity    = 0.125,   # Thermal conductivity [W/m/K]
    density         = 2170.0,  # Density [kg/m³]
    heat_capacity   = 600.0,   # Specific heat capacity [J/kg/K]
    reflectance_vis = 0.059,   # Reflectance in visible light [-]
    reflectance_ir  = 0.0,     # Reflectance in thermal infrared [-]
    emissivity      = 0.9,     # Thermal emissivity [-]
)

grid_params = GridParams(;
    z_max   = 0.6,  # Depth of the lower boundary [m]
    n_depth = 61,   # Number of depth nodes
)

##= Define thermophysical problem =##
problem = BinaryAsteroidThermoPhysicalProblem(
    (shape1, shape2),
    thermo_params,
    grid_params;
    with_self_shadowing      = true,
    with_self_heating        = true,
    with_mutual_shadowing    = true,
    with_mutual_heating      = true,
    upper_boundary_condition = RadiationBoundaryCondition(),
    lower_boundary_condition = InsulationBoundaryCondition(),
)

##= Output specification =##
output_times        = collect(et_range)[end-n_step_in_cycle:end]
subsurface_face_ids = [1, 2, 3, 4, 10]

output = BinaryAsteroidOutputSpec(
    SingleAsteroidOutputSpec(output_times, subsurface_face_ids),
    SingleAsteroidOutputSpec(output_times, subsurface_face_ids),
)

##= Run TPM =##
solution = solve(problem, ExplicitEuler();
    ephem                         = ephem,
    output                        = output,
    initial_temperature_primary   = 200.0,
    initial_temperature_secondary = 200.0,
)

export_solution("path/to/save", solution)
```

## Analyzing Results

After running the TPM, you can access the recorded physical quantities through the `solution` object.

```julia
# Time steps at which quantities were recorded [s]
solution.times

# Power balance [W]
solution.absorbed_power  # Solar power absorbed by the whole surface
solution.emitted_power   # Thermal power emitted from the whole surface

# Surface temperature [K] — matrix of size (n_face, n_time)
solution.surface_temperature

# Subsurface temperature [K] — Dict with face ID as key, matrix (n_depth, n_time) as value
solution.subsurface_temperature[2]  # Subsurface profile at face 2

# For a binary asteroid, results for each body are stored separately:
solution.primary.surface_temperature
solution.secondary.surface_temperature
```

## Computing Thermal Forces and Torques

To compute thermal recoil forces and torques in the inertial frame, provide rotation matrices from the body-fixed frame to the inertial frame when constructing ephemerides.

```julia
using AsteroidShapeModels
using AsteroidThermoPhysicalModels
using Rotations
using SPICE

# ...load kernels, shape model, ephem times as above...

r_sun              = [SPICE.spkpos("SUN", et, "RYUGU_FIXED", "None", "RYUGU")[1] for et in et_range]
R_body_to_inertial = [RotMatrix{3}(SPICE.pxform("RYUGU_FIXED", "J2000", et)) for et in et_range]

r_sun .*= 1000  # Convert [km] to [m]

# Providing R_body_to_inertial enables inertial-frame force/torque output
ephem = SingleAsteroidEphemerides(et_range, r_sun, R_body_to_inertial)

# ...define problem as above...

# Enable force and torque output
output = SingleAsteroidOutputSpec(output_times, subsurface_face_ids;
    save_forces  = true,
    save_torques = true,
)

solution = solve(problem, ExplicitEuler();
    ephem               = ephem,
    output              = output,
    initial_temperature = 200.0,
)

# Net thermal force in the inertial frame [N] — Vector of length n_time
solution.forces

# Net thermal torque in the inertial frame [N⋅m] — Vector of length n_time
solution.torques
```

If only per-face forces in the body-fixed frame are needed (without rotation matrices), use `save_face_forces=true` with a `SingleAsteroidEphemerides{Nothing}` (i.e., without `R_body_to_inertial`):

```julia
ephem  = SingleAsteroidEphemerides(et_range, r_sun)  # no rotation matrix required
output = SingleAsteroidOutputSpec(output_times, subsurface_face_ids; save_face_forces=true)

solution = solve(problem, ExplicitEuler(); ephem=ephem, output=output, initial_temperature=200.0)

# Per-face thermal force in the body-fixed frame [N] — matrix of size (n_face, n_time)
solution.face_forces
```
