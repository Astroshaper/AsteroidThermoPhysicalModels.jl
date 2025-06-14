# Usage Examples

This section provides detailed examples of how to use `AsteroidThermoPhysicalModels.jl` for various scenarios.
For more detailed examples, please refer to the [Astroshaper-examples](https://github.com/Astroshaper/Astroshaper-examples) repository.


## Single Asteroid Example (Ryugu)

This example demonstrates how to set up and run a thermophysical model for asteroid Ryugu using SPICE kernels for ephemerides.

```julia
using AsteroidShapeModels
using AsteroidThermoPhysicalModels
using Downloads
using LinearAlgebra
using Rotations
using SPICE
using StaticArrays

##= Download Files =##
paths_kernel = [
    "lsk/naif0012.tls",
    "pck/hyb2_ryugu_shape_v20190328.tpc",
    "fk/hyb2_ryugu_v01.tf",
    "spk/2162173_Ryugu.bsp",
]
paths_shape = [
    "SHAPE_SFM_49k_v20180804.obj",
]

for path_kernel in paths_kernel
    url_kernel = "https://data.darts.isas.jaxa.jp/pub/hayabusa2/old/2020/spice_bundle/spice_kernels/$(path_kernel)"
    filepath = joinpath("kernel", path_kernel)
    mkpath(dirname(filepath))
    isfile(filepath) || Downloads.download(url_kernel, filepath)
end

for path_shape in paths_shape
    url_shape = "https://data.darts.isas.jaxa.jp/pub/hayabusa2/paper/Watanabe_2019/$(path_shape)"
    filepath = joinpath("shape", path_shape)
    mkpath(dirname(filepath))
    isfile(filepath) || Downloads.download(url_shape, filepath)
end

##= Load data with SPICE =##
for path_kernel in paths_kernel
    filepath = joinpath("kernel", path_kernel)
    SPICE.furnsh(filepath)
end

##= Ephemerides =##
P = SPICE.convrt(7.63262, "hours", "seconds")  # Rotation period of Ryugu

n_cycle = 2  # Number of cycles to perform TPM
n_step_in_cycle = 72  # Number of time steps in one rotation period

et_begin = SPICE.utc2et("2018-07-01T00:00:00")  # Start time of TPM
et_end   = et_begin + P * n_cycle  # End time of TPM
et_range = range(et_begin, et_end; length=n_step_in_cycle*n_cycle+1)

"""
- `time` : Ephemeris times
- `sun`  : Sun's position in the RYUGU_FIXED frame
"""
ephem = (
    time = collect(et_range),
    sun  = [SVector{3}(SPICE.spkpos("SUN", et, "RYUGU_FIXED", "None", "RYUGU")[1]) * 1000 for et in et_range],
)

SPICE.kclear()

##= Load obj file =##
path_obj = joinpath("shape", "SHAPE_SFM_49k_v20180804.obj")
    
shape = load_shape_obj(path_obj; scale=1000, find_visible_facets=true)
n_face = length(shape.faces)  # Number of faces

##= Thermal properties =##
k  = 0.1     # Thermal conductivity [W/m/K]
ρ  = 1270.0  # Density [kg/m³]
Cₚ = 600.0   # Heat capacity [J/kg/K]
    
l = AsteroidThermoPhysicalModels.thermal_skin_depth(P, k, ρ, Cₚ)  # Thermal skin depth [m]
Γ = AsteroidThermoPhysicalModels.thermal_inertia(k, ρ, Cₚ)        # Thermal inertia [tiu]

R_vis = 0.04  # Reflectance in visible light [-]
R_ir  = 0.0   # Reflectance in thermal infrared [-]
ε     = 1.0   # Emissivity [-]

z_max = 0.6   # Depth of the lower boundary of a heat conduction equation [m]
n_depth = 41  # Number of depth steps
Δz = z_max / (n_depth - 1)  # Depth step width [m]

thermo_params = AsteroidThermoPhysicalModels.ThermoParams(P, l, Γ, R_vis, R_ir, ε, z_max, Δz, n_depth)

##= Setting of TPM =##
stpm = AsteroidThermoPhysicalModels.SingleAsteroidTPM(shape, thermo_params;
    SELF_SHADOWING = true,
    SELF_HEATING   = true,
    SOLVER         = AsteroidThermoPhysicalModels.CrankNicolsonSolver(thermo_params),
    BC_UPPER       = AsteroidThermoPhysicalModels.RadiationBoundaryCondition(),
    BC_LOWER       = AsteroidThermoPhysicalModels.InsulationBoundaryCondition(),
)
AsteroidThermoPhysicalModels.init_temperature!(stpm, 200)

##= Run TPM =##
times_to_save = ephem.time[end-n_step_in_cycle:end]  # Save temperature during the final rotation
face_ID = [1, 2, 3, 4, 10]  # Face indices to save subsurface temperature

result = AsteroidThermoPhysicalModels.run_TPM!(stpm, ephem, times_to_save, face_ID)
AsteroidThermoPhysicalModels.export_TPM_results("path/to/save", result)
```

## Binary Asteroid Example (Didymos-Dimorphos)

This example demonstrates how to set up and run a thermophysical model for the binary asteroid system Didymos-Dimorphos using SPICE kernels for ephemerides.

```julia
using AsteroidShapeModels
using AsteroidThermoPhysicalModels
using Downloads
using LinearAlgebra
using Rotations
using SPICE
using StaticArrays

##= SPICE kernels =##
paths_kernel = [
    "fk/hera_v10.tf",
    "lsk/naif0012.tls",
    "pck/hera_didymos_v06.tpc",
    "spk/de432s.bsp",
    "spk/didymos_hor_000101_500101_v01.bsp",
    "spk/didymos_gmv_260901_311001_v01.bsp",
]

##= Shape models =##
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
    url_kernel = "https://s2e2.cosmos.esa.int/bitbucket/projects/SPICE_KERNELS/repos/hera/raw/kernels/dsk/$(path_shape)?at=refs%2Ftags%2Fv161_20230929_001"
    filepath = joinpath("shape", path_shape)
    mkpath(dirname(filepath))
    isfile(filepath) || Downloads.download(url_kernel, filepath)
end

##= Load the SPICE kernels =##
for path_kernel in paths_kernel
    filepath = joinpath("kernel", path_kernel)
    SPICE.furnsh(filepath)
end

##= Ephemerides =##
P₁ = SPICE.convrt(2.2593, "hours", "seconds")  # Rotation period of Didymos
P₂ = SPICE.convrt(11.93 , "hours", "seconds")  # Rotation period of Dimorphos

n_cycle = 2  # Number of cycles to perform TPM
n_step_in_cycle = 72  # Number of time steps in one rotation period

et_begin = SPICE.utc2et("2027-02-18T00:00:00")  # Start time of TPM
et_end   = et_begin + P₂ * n_cycle  # End time of TPM
et_range = range(et_begin, et_end; length=n_step_in_cycle*n_cycle+1)

"""
- `time` : Ephemeris times
- `sun1` : Sun's position in the primary's frame (DIDYMOS_FIXED)
- `sun2` : Sun's position in the secondary's frame (DIMORPHOS_FIXED)
- `sec`  : Secondary's position in the primary's frame (DIDYMOS_FIXED)
- `P2S`  : Rotation matrix from primary to secondary frames
- `S2P`  : Rotation matrix from secondary to primary frames
"""
ephem = (
    time = collect(et_range),
    sun1 = [SVector{3}(SPICE.spkpos("SUN"      , et, "DIDYMOS_FIXED"  , "None", "DIDYMOS"  )[1]) * 1000 for et in et_range],
    sun2 = [SVector{3}(SPICE.spkpos("SUN"      , et, "DIMORPHOS_FIXED", "None", "DIMORPHOS")[1]) * 1000 for et in et_range],
    sec  = [SVector{3}(SPICE.spkpos("DIMORPHOS", et, "DIDYMOS_FIXED"  , "None", "DIDYMOS"  )[1]) * 1000 for et in et_range],
    P2S  = [RotMatrix{3}(SPICE.pxform("DIDYMOS_FIXED"  , "DIMORPHOS_FIXED", et)) for et in et_range],
    S2P  = [RotMatrix{3}(SPICE.pxform("DIMORPHOS_FIXED", "DIDYMOS_FIXED"  , et)) for et in et_range],
)

SPICE.kclear()

##= Load the shape models =##
path_shape1_obj = joinpath("shape", "g_50677mm_rad_obj_didy_0000n00000_v001.obj")
path_shape2_obj = joinpath("shape", "g_08438mm_lgt_obj_dimo_0000n00000_v002.obj")
    
shape1 = load_shape_obj(path_shape1_obj; scale=1000, find_visible_facets=true)
shape2 = load_shape_obj(path_shape2_obj; scale=1000, find_visible_facets=true)

n_face_shape1 = length(shape1.faces)  # Number of faces of Didymos
n_face_shape2 = length(shape2.faces)  # Number of faces of Dimorphos
    
##= Thermal properties of Didymos & Dimorphos [cf. Michel+2016; Naidu+2020] =##
k  = 0.125   # Thermal conductivity [W/m/K]
ρ  = 2170.0  # Density [kg/m³]
Cₚ = 600.0   # Heat capacity [J/kg/K]

l₁ = AsteroidThermoPhysicalModels.thermal_skin_depth(P₁, k, ρ, Cₚ)  # Thermal skin depth for Didymos
l₂ = AsteroidThermoPhysicalModels.thermal_skin_depth(P₂, k, ρ, Cₚ)  # Thermal skin depth for Dimorphos
Γ = AsteroidThermoPhysicalModels.thermal_inertia(k, ρ, Cₚ)          # Thermal inertia for Didymos and Dimorphos [tiu]

R_vis = 0.059  # Reflectance in visible light [-]
R_ir  = 0.0    # Reflectance in thermal infrared [-]
ε     = 0.9    # Emissivity [-]

z_max = 0.6   # Depth of the lower boundary of a heat conduction equation [m]
n_depth = 41  # Number of depth steps
Δz = z_max / (n_depth - 1)  # Depth step width [m]

thermo_params1 = AsteroidThermoPhysicalModels.ThermoParams(P₁, l₁, Γ, R_vis, R_ir, ε, z_max, Δz, n_depth)
thermo_params2 = AsteroidThermoPhysicalModels.ThermoParams(P₂, l₂, Γ, R_vis, R_ir, ε, z_max, Δz, n_depth)

##= Setting of TPM =##
stpm1 = AsteroidThermoPhysicalModels.SingleAsteroidTPM(shape1, thermo_params1;
    SELF_SHADOWING = true,
    SELF_HEATING   = true,
    SOLVER         = AsteroidThermoPhysicalModels.CrankNicolsonSolver(thermo_params1),
    BC_UPPER       = AsteroidThermoPhysicalModels.RadiationBoundaryCondition(),
    BC_LOWER       = AsteroidThermoPhysicalModels.InsulationBoundaryCondition(),
)

stpm2 = AsteroidThermoPhysicalModels.SingleAsteroidTPM(shape2, thermo_params2;
    SELF_SHADOWING = true,
    SELF_HEATING   = true,
    SOLVER         = AsteroidThermoPhysicalModels.CrankNicolsonSolver(thermo_params2),
    BC_UPPER       = AsteroidThermoPhysicalModels.RadiationBoundaryCondition(),
    BC_LOWER       = AsteroidThermoPhysicalModels.InsulationBoundaryCondition(),
)

btpm  = AsteroidThermoPhysicalModels.BinaryAsteroidTPM(stpm1, stpm2; MUTUAL_SHADOWING=true, MUTUAL_HEATING=true)
AsteroidThermoPhysicalModels.init_temperature!(btpm, 200.)
    
##= Run TPM =##
times_to_save = ephem.time[end-n_step_in_cycle:end]  # Save temperature during the final rotation
face_ID_pri = [1, 2, 3, 4, 10]  # Face indices to save subsurface temperature of the primary
face_ID_sec = [1, 2, 3, 4, 20]  # Face indices to save subsurface temperature of the secondary

result = AsteroidThermoPhysicalModels.run_TPM!(btpm, ephem, times_to_save, face_ID_pri, face_ID_sec)
AsteroidThermoPhysicalModels.export_TPM_results("path/to/save", result)
```

## Analyzing Results

After running the TPM, you can analyze the results to understand the thermal behavior of the asteroid and its non-gravitational effects.

```julia
# Access physical quantities
result.times  # Time steps for thermophysical modeling [s]
result.E_in   # Input energy per second on the whole surface [W]
result.E_out  # Output energy per second from the whole surface [W]
result.force  # Thermal force on the asteroid at the body-fixed frame [N]
result.torque # Thermal torque on the asteroid at the body-fixed frame [N ⋅ m]

# For a binary asteroid, you can access the results of the primary and secondary as follows:
result.pri.force  # Thermal force on the primary body at the body-fixed frame [N]
result.sec.force  # Thermal force on the secondary body at the body-fixed frame [N]

# Access surface temperature [K]
# A matrix in size of `(n_face, n_time)`
# - `n_face` : Number of faces of the shape model
# - `n_time` : Number of time steps to save surface temperature
result.surface_temperature

# Access subsurface temperatures [K]
# As a function of depth [m] and time [s], `Dict` with face ID as key and a matrix `(n_depth, n_time)` as an entry.
# - `n_depth` : The number of the depth nodes
# - `n_time`  : The number of time steps to save temperature
result.subsurface_temperature[2]  # Subsurface temperature for face ID 2
```
