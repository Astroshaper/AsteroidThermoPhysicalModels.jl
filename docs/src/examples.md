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

output = SingleAsteroidOutputSpec(output_times; subsurface_face_ids)

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
    SingleAsteroidOutputSpec(output_times; subsurface_face_ids),
    SingleAsteroidOutputSpec(output_times; subsurface_face_ids),
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
output = SingleAsteroidOutputSpec(output_times;
    subsurface_face_ids = subsurface_face_ids,
    save_forces         = true,
    save_torques        = true,
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
output = SingleAsteroidOutputSpec(output_times; subsurface_face_ids, save_face_forces=true)

solution = solve(problem, ExplicitEuler(); ephem=ephem, output=output, initial_temperature=200.0)

# Per-face thermal force in the body-fixed frame [N] — matrix of size (n_face, n_time)
solution.face_forces
```

## Surface Roughness

Surface roughness is modelled by attaching a small shape model — a *roughness model*, for example a spherical crater — to the facets of a `ShapeModel` from `AsteroidShapeModels.jl` (stored in its `roughness` field). Each facet with a roughness model is then solved as a full thermophysical model of its own, in the local frame of the facet, with self-shadowing and self-heating inside the roughness model. This is the origin of thermal-infrared beaming.

```julia
using AsteroidShapeModels
using AsteroidThermoPhysicalModels

# Load the global shape
shape = load_shape_obj("path/to/shape.obj"; scale=1000)

# A spherical crater of radius 0.4 and depth 0.1 (in the units of the roughness model),
# discretised on an 8 × 8 grid, attached to every facet
crater = create_shape_crater(0.4, 0.1; Nx=8, Ny=8)
add_roughness_models!(shape, crater)
# ...or to selected facets only: add_roughness_models!(shape, crater, face_idx)

# The problem, ephemerides, output and solver are defined exactly as for a plain ShapeModel
problem = SingleAsteroidThermoPhysicalProblem(shape, thermo_params, grid_params;
    with_self_shadowing = true,
    with_self_heating   = true,
)
solution = solve(problem, CrankNicolson(); ephem=ephem, output=output, initial_temperature=200.0)
```

The recorded `surface_temperature` and `subsurface_temperature` are those of the global facets, which are solved independently of their roughness models and serve as the smooth-surface baseline of the same run. `face_forces`, `forces` and `torques` include the roughness: on a facet with a roughness model the force is the sum over its sub-facets, counted for the area of the facet (see *Facets with a roughness model* in the physical model). `absorbed_power` and `emitted_power` count such facets from their sub-facets in the same way.

The problem's `with_self_heating` governs the global facets and, through them, the irradiation the sub-facets receive from the rest of the body; self-shadowing and self-heating inside a roughness model are always on.

To record the temperatures of the rough surface itself, list the facets whose roughness models to save:

```julia
output = SingleAsteroidOutputSpec(output_times;
    subsurface_face_ids = subsurface_face_ids,
    roughness_face_ids  = [1, 7],
)
solution = solve(problem, CrankNicolson(); ephem=ephem, output=output, initial_temperature=200.0)

solution.roughness_surface_temperature[7]   # (n_sub, n_output): sub-facet j of the crater on facet 7, at each output time
```

`export_solution` then also writes `roughness_surface_temperature.csv` in long format (`time`, `face_id`, `sub_face_id`, `temperature`), since roughness models may differ in size from facet to facet.

From these temperatures, the radiance or brightness temperature of every facet towards an observer is a post-processing step, so any number of observer directions can be tried on one solution. For a thermal-infrared image, list every facet in `roughness_face_ids`:

```julia
d̂_obs = normalize(r_observer)   # direction to the observer in the body-fixed frame
i_save = 3                        # index into output_times, e.g. the image epoch

L   = directional_radiance(problem, solution, i_save, d̂_obs)            # W/m²/sr, one value per facet
T_b = brightness_temperature(problem, solution, i_save, d̂_obs)          # K
L_λ = directional_radiance(problem, solution, i_save, d̂_obs; λ=10e-6)   # spectral, at 10 μm
```

Facets with a recorded roughness model radiate according to their sub-facet temperatures and the crater walls in view (see *Direction-Dependent Radiance of a Rough Facet* in the physical model); the others as smooth Lambertian surfaces. Facets seen from behind give `NaN`. The per-facet vector is ready to be placed on an image by a ray-caster such as `FOVSimulator.jl`.
