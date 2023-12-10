var documenterSearchIndex = {"docs":
[{"location":"api/#API","page":"API","title":"API","text":"","category":"section"},{"location":"api/","page":"API","title":"API","text":"Modules = [AsteroidThermoPhysicalModels]","category":"page"},{"location":"api/#AsteroidThermoPhysicalModels.CERES","page":"API","title":"AsteroidThermoPhysicalModels.CERES","text":"Ceres\n\n\n\n\n\n","category":"constant"},{"location":"api/#AsteroidThermoPhysicalModels.DIDYMOS","page":"API","title":"AsteroidThermoPhysicalModels.DIDYMOS","text":"Asteroid 65803 Didymos (1996 GT)\n\nPhysical parameters\n\n:GM : GM\n:M  : Mass\n:μ  : Standard gravitational parameter about Sun and Ryugu\n\nOrbital elements\n\n:a  : Semi-mojor axis [AU]\n:e  : Eccentricity [-]\n:I  : Inclination [deg]\n:Ω  : Longitude of the ascending node [deg]\n:ω  : Argument of periapsis [deg]\n:Φ  : # Mean anomaly [deg]\n\nSpin parameters\n\n:α  : Right ascension (RA) in equatorial coordinate system [deg]\n:δ  : Declination (Dec) in equatorial coordinate system [deg]\n:P  : Rotation period [h]\n\nReferences\n\nSmall Body Database Lookup: https://ssd.jpl.nasa.gov/tools/sbdb_lookup.html#/?sstr=didymos\nNaidu et al. (2020)\n\n\n\n\n\n","category":"constant"},{"location":"api/#AsteroidThermoPhysicalModels.DIMORPHOS","page":"API","title":"AsteroidThermoPhysicalModels.DIMORPHOS","text":"Dimorphos\n\n:a : Semi-major axis of the mutual orbit with Dimorphos [m]\n:P : Rotation periof of the mutual orbit with Dimorphos [h]\n\n\n\n\n\n","category":"constant"},{"location":"api/#AsteroidThermoPhysicalModels.EARTH","page":"API","title":"AsteroidThermoPhysicalModels.EARTH","text":"Earth\n\n\n\n\n\n","category":"constant"},{"location":"api/#AsteroidThermoPhysicalModels.ERIS","page":"API","title":"AsteroidThermoPhysicalModels.ERIS","text":"Eris\n\n\n\n\n\n","category":"constant"},{"location":"api/#AsteroidThermoPhysicalModels.JUPITER","page":"API","title":"AsteroidThermoPhysicalModels.JUPITER","text":"Jupiter\n\n\n\n\n\n","category":"constant"},{"location":"api/#AsteroidThermoPhysicalModels.MARS","page":"API","title":"AsteroidThermoPhysicalModels.MARS","text":"Mars\n\n\n\n\n\n","category":"constant"},{"location":"api/#AsteroidThermoPhysicalModels.MERCURY","page":"API","title":"AsteroidThermoPhysicalModels.MERCURY","text":"Mercury\n\n\n\n\n\n","category":"constant"},{"location":"api/#AsteroidThermoPhysicalModels.MOON","page":"API","title":"AsteroidThermoPhysicalModels.MOON","text":"Moon\n\n\n\n\n\n","category":"constant"},{"location":"api/#AsteroidThermoPhysicalModels.NEPTUNE","page":"API","title":"AsteroidThermoPhysicalModels.NEPTUNE","text":"Neptune\n\n\n\n\n\n","category":"constant"},{"location":"api/#AsteroidThermoPhysicalModels.PLUTO","page":"API","title":"AsteroidThermoPhysicalModels.PLUTO","text":"Pluto\n\n\n\n\n\n","category":"constant"},{"location":"api/#AsteroidThermoPhysicalModels.RYUGU","page":"API","title":"AsteroidThermoPhysicalModels.RYUGU","text":"Asteroid 162173 Ryugu\n\nPhysical parameters\n\n:GM : GM\n:M  : Mass\n:μ  : Standard gravitational parameter about Sun and Ryugu\n\nOrbital elements\n\n:a  : Semi-mojor axis [AU]\n:e  : Eccentricity [-]\n:I  : Inclination [deg]\n:Ω  : Longitude of the ascending node [deg]\n:ω  : Argument of periapsis [deg]\n:Φ  : # Mean anomaly [deg]\n\nSpin parameters\n\n:α  : Right ascension (RA) in equatorial coordinate system [deg]\n:δ  : Declination (Dec) in equatorial coordinate system [deg]\n:P  : Rotation period [h]\n\n\n\n\n\n","category":"constant"},{"location":"api/#AsteroidThermoPhysicalModels.SATURN","page":"API","title":"AsteroidThermoPhysicalModels.SATURN","text":"Saturn\n\n\n\n\n\n","category":"constant"},{"location":"api/#AsteroidThermoPhysicalModels.URANUS","page":"API","title":"AsteroidThermoPhysicalModels.URANUS","text":"Uranus\n\n\n\n\n\n","category":"constant"},{"location":"api/#AsteroidThermoPhysicalModels.VENUS","page":"API","title":"AsteroidThermoPhysicalModels.VENUS","text":"Venus\n\n\n\n\n\n","category":"constant"},{"location":"api/#AsteroidThermoPhysicalModels.BackwardEulerSolver","page":"API","title":"AsteroidThermoPhysicalModels.BackwardEulerSolver","text":"Type of the backward Euler method:\n\nImplicit in time (Unconditionally stable in the heat conduction equation)\nFirst order in time\n\nThe BackwardEulerSolver type has vectors for the tridiagonal matrix algorithm.\n\n\n\n\n\n","category":"type"},{"location":"api/#AsteroidThermoPhysicalModels.BinaryTPM","page":"API","title":"AsteroidThermoPhysicalModels.BinaryTPM","text":"struct BinaryTPM <: ThermoPhysicalModel\n\nFields\n\npri              : TPM for the primary\nsec              : TPM for the secondary\nMUTUAL_SHADOWING : Flag to consider mutual shadowing\nMUTUAL_HEATING   : Flag to consider mutual heating\n\n\n\n\n\n","category":"type"},{"location":"api/#AsteroidThermoPhysicalModels.BinaryTPM-Tuple{Any, Any}","page":"API","title":"AsteroidThermoPhysicalModels.BinaryTPM","text":"BinaryTPM(pri, sec; MUTUAL_SHADOWING=true, MUTUAL_HEATING=true) -> btpm\n\nConstruct a thermophysical model for a binary asteroid (BinaryTPM).\n\n\n\n\n\n","category":"method"},{"location":"api/#AsteroidThermoPhysicalModels.BinaryTPMResult","page":"API","title":"AsteroidThermoPhysicalModels.BinaryTPMResult","text":"struct BinaryTPMResult\n\nOutput data format for BinaryTPM\n\nFields\n\npri : TPM result for the primary\nsec : TPM result for the secondary\n\n\n\n\n\n","category":"type"},{"location":"api/#AsteroidThermoPhysicalModels.BinaryTPMResult-Tuple{AsteroidThermoPhysicalModels.BinaryTPM, Any, Vector{Float64}, Vector{Int64}, Vector{Int64}}","page":"API","title":"AsteroidThermoPhysicalModels.BinaryTPMResult","text":"Outer constructor of BinaryTPMResult\n\nArguments\n\nbtpm          : Thermophysical model for a binary asteroid\nephem         : Ephemerides\ntimes_to_save : Timesteps to save temperature (Common to both the primary and the secondary)\nface_ID_pri   : Face indices to save subsurface temperature of the primary\nface_ID_sec   : Face indices to save subsurface temperature of the secondary\n\n\n\n\n\n","category":"method"},{"location":"api/#AsteroidThermoPhysicalModels.BoundaryCondition","page":"API","title":"AsteroidThermoPhysicalModels.BoundaryCondition","text":"Abstract type of a boundary condition for a heat conduction equation\n\n\n\n\n\n","category":"type"},{"location":"api/#AsteroidThermoPhysicalModels.CrankNicolsonSolver","page":"API","title":"AsteroidThermoPhysicalModels.CrankNicolsonSolver","text":"Type of the Crank-Nicolson method:\n\nImplicit in time (Unconditionally stable in the heat conduction equation)\nSecond order in time\n\nThe CrankNicolsonSolver type has vectors for the tridiagonal matrix algorithm.\n\nReferences\n\nhttps://en.wikipedia.org/wiki/Crank–Nicolson_method\n\n\n\n\n\n","category":"type"},{"location":"api/#AsteroidThermoPhysicalModels.ForwardEulerSolver","page":"API","title":"AsteroidThermoPhysicalModels.ForwardEulerSolver","text":"Type of the forward Euler method:\n\nExplicit in time\nFirst order in time\n\nThe ForwardEulerSolver type includes a vector for the temperature at the next time step.\n\n\n\n\n\n","category":"type"},{"location":"api/#AsteroidThermoPhysicalModels.HeatConductionSolver","page":"API","title":"AsteroidThermoPhysicalModels.HeatConductionSolver","text":"Abstract type of a solver for a heat conduction equation \n\n\n\n\n\n","category":"type"},{"location":"api/#AsteroidThermoPhysicalModels.InsulationBoundaryCondition","page":"API","title":"AsteroidThermoPhysicalModels.InsulationBoundaryCondition","text":"Singleton type of insulation boundary condition\n\n\n\n\n\n","category":"type"},{"location":"api/#AsteroidThermoPhysicalModels.IsothermalBoundaryCondition","page":"API","title":"AsteroidThermoPhysicalModels.IsothermalBoundaryCondition","text":"Type of isothermal boundary condition\n\n\n\n\n\n","category":"type"},{"location":"api/#AsteroidThermoPhysicalModels.NonUniformThermoParams","page":"API","title":"AsteroidThermoPhysicalModels.NonUniformThermoParams","text":"struct NonUniformThermoParams\n\nFields\n\nP     : Cycle of thermal cycle (rotation period) [sec]\nl     : Thermal skin depth [m]\nΓ     : Thermal inertia [J ⋅ m⁻² ⋅ K⁻¹ ⋅ s⁻⁰⁵ (tiu)]\nA_B   : Bond albedo\nA_TH  : Albedo at thermal radiation wavelength\nε     : Emissivity\nz_max : Depth of the bottom of a heat conduction equation [m]\nΔz    : Depth step width [m]\nNz    : Number of depth steps\n\n\n\n\n\n","category":"type"},{"location":"api/#AsteroidThermoPhysicalModels.RadiationBoundaryCondition","page":"API","title":"AsteroidThermoPhysicalModels.RadiationBoundaryCondition","text":"Singleton type of radiation boundary condition\n\n\n\n\n\n","category":"type"},{"location":"api/#AsteroidThermoPhysicalModels.ShapeModel","page":"API","title":"AsteroidThermoPhysicalModels.ShapeModel","text":"ShapeModel\n\nA polyhedral shape model of an asteroid.\n\nFields\n\nnodes         : Vector of node positions\nfaces         : Vector of vertex indices of faces\nface_centers  : Center position of each face\nface_normals  : Normal vector of each face\nface_areas    : Area of of each face\nvisiblefacets : Vector of vector of VisibleFacet\n\n\n\n\n\n","category":"type"},{"location":"api/#AsteroidThermoPhysicalModels.SingleTPM","page":"API","title":"AsteroidThermoPhysicalModels.SingleTPM","text":"struct SingleTPM <: ThermoPhysicalModel\n\nFields\n\nshape          : Shape model\nthermo_params  : Thermophysical parameters\nflux           : Flux on each face. Matrix of size (Number of faces, 3). Three components are:\nflux[:, 1]     : F_sun,  flux of direct sunlight\nflux[:, 2]     : F_scat, flux of scattered light\nflux[:, 3]     : F_rad,  flux of thermal emission from surrounding surface\ntemperature    : Temperature matrix (Nz, Ns) according to the number of depth cells Nz and the number of faces Ns.\nface_forces    : Thermal force on each face\nforce          : Thermal recoil force at body-fixed frame (Yarkovsky effect)\ntorque         : Thermal recoil torque at body-fixed frame (YORP effect)\nSELF_SHADOWING : Flag to consider self-shadowing\nSELF_HEATING   : Flag to consider self-heating\nSOLVER         : Solver of heat conduction equation\nBC_UPPER       : Boundary condition at the upper boundary\nBC_LOWER       : Boundary condition at the lower boundary\n\nTO DO:\n\nroughness_maps   ::ShapeModel[]\n\n\n\n\n\n","category":"type"},{"location":"api/#AsteroidThermoPhysicalModels.SingleTPM-Tuple{Any, Any}","page":"API","title":"AsteroidThermoPhysicalModels.SingleTPM","text":"SingleTPM(shape, thermo_params; SELF_SHADOWING=true, SELF_HEATING=true) -> stpm\n\nConstruct a thermophysical model for a single asteroid (SingleTPM).\n\nArguments\n\nshape          : Shape model\nthermo_params  : Thermophysical parameters\n\nKeyword arguments\n\nSELF_SHADOWING : Flag to consider self-shadowing\nSELF_HEATING   : Flag to consider self-heating\nSOLVER         : Solver of heat conduction equation\nBC_UPPER       : Boundary condition at the upper boundary\nBC_LOWER       : Boundary condition at the lower boundary\n\n\n\n\n\n","category":"method"},{"location":"api/#AsteroidThermoPhysicalModels.SingleTPMResult","page":"API","title":"AsteroidThermoPhysicalModels.SingleTPMResult","text":"struct SingleTPMResult\n\nOutput data format for SingleTPM \n\nFields\n\nSaved at all time steps\n\ntimes  : Timesteps, given the same vector as ephem.time [s]\nE_in   : Input energy per second on the whole surface [W]\nE_out  : Output enegey per second from the whole surface [W]\nE_cons : Energy conservation ratio [-], ratio of total energy going out to total energy coming in in the last rotation cycle\nforce  : Thermal force on the asteroid [N]\ntorque : Thermal torque on the asteroid [N ⋅ m]\n\nSaved only at the time steps desired by the user\n\ntimes_to_save : Timesteps to save temperature [s]\ndepth_nodes   : Depths of the calculation nodes for 1-D heat conduction [m], a vector of size Nz\nsurface_temperature     : Surface temperature [K], a matrix in size of (Ns, Nt).\nNs : Number of faces\nNt : Number of time steps to save surface temperature\nsubsurface_temperature     : Temperature [K] as a function of depth [m] and time [s], Dict with face ID as key and a matrix (Nz, Nt) as an entry.\nNz : The number of the depth nodes\nNt : The number of time steps to save temperature\n\n\n\n\n\n","category":"type"},{"location":"api/#AsteroidThermoPhysicalModels.SingleTPMResult-Tuple{AsteroidThermoPhysicalModels.SingleTPM, Any, Vector{Float64}, Vector{Int64}}","page":"API","title":"AsteroidThermoPhysicalModels.SingleTPMResult","text":"Outer constructor of SingleTPMResult\n\nArguments\n\nstpm          : Thermophysical model for a single asteroid\nephem         : Ephemerides\ntimes_to_save : Timesteps to save temperature\nface_ID       : Face indices to save subsurface temperature\n\n\n\n\n\n","category":"method"},{"location":"api/#AsteroidThermoPhysicalModels.ThermoPhysicalModel","page":"API","title":"AsteroidThermoPhysicalModels.ThermoPhysicalModel","text":"Abstract type of a thermophysical model\n\n\n\n\n\n","category":"type"},{"location":"api/#AsteroidThermoPhysicalModels.UniformThermoParams","page":"API","title":"AsteroidThermoPhysicalModels.UniformThermoParams","text":"struct UniformThermoParams\n\nFields\n\nP     : Thermal cycle (rotation period) [sec]\nl     : Thermal skin depth [m]\nΓ     : Thermal inertia [J ⋅ m⁻² ⋅ K⁻¹ ⋅ s⁻⁰⁵ (tiu)]\nA_B   : Bond albedo\nA_TH  : Albedo at thermal radiation wavelength\nε     : Emissivity\nz_max : Depth of the bottom of a heat conduction equation [m]\nΔz    : Depth step width [m]\nNz    : Number of depth steps\n\n\n\n\n\n","category":"type"},{"location":"api/#AsteroidThermoPhysicalModels.VisibleFacet","page":"API","title":"AsteroidThermoPhysicalModels.VisibleFacet","text":"struct VisibleFacet\n\nIndex of an interfacing facet and its view factor\n\nFields\n\nid : Index of the interfacing facet\nf  : View factor from facet i to j\nd  : Distance from facet i to j\nd̂  : Normal vector from facet i to j\n\n\n\n\n\n","category":"type"},{"location":"api/#AsteroidThermoPhysicalModels.backward_euler!-Tuple{AsteroidThermoPhysicalModels.SingleTPM, Any}","page":"API","title":"AsteroidThermoPhysicalModels.backward_euler!","text":"backward_euler!(stpm::SingleTPM, Δt)\n\nPredict the temperature at the next time step by the backward Euler method.\n\nImplicit in time (Unconditionally stable in the heat conduction equation)\nFirst order in time\nSecond order in space\n\nIn this function, the heat conduction equation is non-dimensionalized in time and length.\n\n\n\n\n\n","category":"method"},{"location":"api/#AsteroidThermoPhysicalModels.concave_spherical_segment-NTuple{6, Real}","page":"API","title":"AsteroidThermoPhysicalModels.concave_spherical_segment","text":"concave_spherical_segment(r, h, xc, yc, x, y) -> z\n\nReturn the z-coordinate of a concave spherical segment.\n\nArguments\n\nr : Crater radius\nh : Crater depth\nxc: x-coordinate of crater center\nyc: y-coordinate of crater center\nx : x-coordinate where to calculate z\ny : y-coordinate where to calculate z\n\n\n\n\n\n","category":"method"},{"location":"api/#AsteroidThermoPhysicalModels.concave_spherical_segment-Tuple{Real, Real}","page":"API","title":"AsteroidThermoPhysicalModels.concave_spherical_segment","text":"concave_spherical_segment(r, h; Nx=2^5, Ny=2^5, xc=0.5, yc=0.5) -> xs, ys, zs\n\nReturn (x, y, z) grid of a concave spherical segment.\n\nArguments\n\nr : Crater radius\nh : Crater depth\nNx: Number of nodes in the x-direction\nNy: Number of nodes in the y-direction\nxc: x-coordinate of crater center\nyc: y-coordinate of crater center\n\n\n\n\n\n","category":"method"},{"location":"api/#AsteroidThermoPhysicalModels.crank_nicolson!-Tuple{AsteroidThermoPhysicalModels.SingleTPM, Any}","page":"API","title":"AsteroidThermoPhysicalModels.crank_nicolson!","text":"crank_nicolson!(stpm::SingleTPM, Δt)\n\nPredict the temperature at the next time step by the Crank-Nicolson method.\n\nImplicit in time (Unconditionally stable in the heat conduction equation)\nSecond order in time\nSecond order in space\n\nIn this function, the heat conduction equation is non-dimensionalized in time and length.\n\n\n\n\n\n","category":"method"},{"location":"api/#AsteroidThermoPhysicalModels.crater_curvature_radius-Tuple{Real, Real}","page":"API","title":"AsteroidThermoPhysicalModels.crater_curvature_radius","text":"crater_curvature_radius(r, h) -> R\n\nReturn the curvature radius of a concave spherical segment.\n\nArguments\n\nr: Crater radius\nh: Crater depth\n\n\n\n\n\n","category":"method"},{"location":"api/#AsteroidThermoPhysicalModels.energy_in-Tuple{AsteroidThermoPhysicalModels.SingleTPM}","page":"API","title":"AsteroidThermoPhysicalModels.energy_in","text":"energy_in(stpm::SingleTPM) -> E_in\n\nInput energy per second on the whole surface [W]\n\n\n\n\n\n","category":"method"},{"location":"api/#AsteroidThermoPhysicalModels.energy_out-Tuple{AsteroidThermoPhysicalModels.SingleTPM}","page":"API","title":"AsteroidThermoPhysicalModels.energy_out","text":"energy_out(stpm::SingleTPM) -> E_out\n\nOutput enegey per second from the whole surface [W]\n\nArguments\n\nstpm : Thermophysical model for a single asteroid\n\n\n\n\n\n","category":"method"},{"location":"api/#AsteroidThermoPhysicalModels.export_TPM_results-Tuple{Any, AsteroidThermoPhysicalModels.BinaryTPMResult}","page":"API","title":"AsteroidThermoPhysicalModels.export_TPM_results","text":"export_TPM_results(filepath, result::BinaryTPMResult)\n\nExport the result of BinaryTPM to CSV files.\n\nArguments\n\ndirpath : Path to the directory to save CSV files\nresult  : Output data format for BinaryTPM\n\n\n\n\n\n","category":"method"},{"location":"api/#AsteroidThermoPhysicalModels.export_TPM_results-Tuple{Any, AsteroidThermoPhysicalModels.SingleTPMResult}","page":"API","title":"AsteroidThermoPhysicalModels.export_TPM_results","text":"export_TPM_results(dirpath, result::SingleTPMResult)\n\nExport the result of SingleTPM to CSV files.\n\nArguments\n\ndirpath :  Path to the directory to save CSV files\nresult  : Output data format for SingleTPM\n\n\n\n\n\n","category":"method"},{"location":"api/#AsteroidThermoPhysicalModels.find_visiblefacets!-Tuple{ShapeModel}","page":"API","title":"AsteroidThermoPhysicalModels.find_visiblefacets!","text":"find_visiblefacets!(obs::Facet, facets)\n\nFind facets that is visible from the facet where the observer is located.\n\nParameters\n\nobs    : Facet where the observer stands\nfacets : Array of Facet\n\n\n\n\n\n","category":"method"},{"location":"api/#AsteroidThermoPhysicalModels.flux_total-NTuple{5, Any}","page":"API","title":"AsteroidThermoPhysicalModels.flux_total","text":"flux_total(A_B, A_TH, F_sun, F_scat, F_rad) -> F_total\n\nTotal energy absorbed by the surface [W/m²]\n\nArguments\n\nA_B    : Bond albedo [-]\nA_TH   : Albedo at thermal infrared wavelength [-]\nF_sun  : Flux of direct sunlight [W/m²]\nF_scat : Flux of scattered light [W/m²]\nF_rad  : Flux of thermal radiation from surrounding surface [W/m²]\n\n\n\n\n\n","category":"method"},{"location":"api/#AsteroidThermoPhysicalModels.forward_euler!-Tuple{AsteroidThermoPhysicalModels.SingleTPM, Any}","page":"API","title":"AsteroidThermoPhysicalModels.forward_euler!","text":"forward_euler!(stpm::SingleTPM, Δt)\n\nPredict the temperature at the next time step by the forward Euler method.\n\nExplicit in time\nFirst order in time\n\nIn this function, the heat conduction equation is non-dimensionalized in time and length.\n\nArguments\n\nstpm : Thermophysical model for a single asteroid\nΔt   : Time step [sec]\n\n\n\n\n\n","category":"method"},{"location":"api/#AsteroidThermoPhysicalModels.grid_to_faces-Tuple{AbstractVector, AbstractVector, AbstractMatrix}","page":"API","title":"AsteroidThermoPhysicalModels.grid_to_faces","text":"grid_to_faces(xs::AbstractVector, ys::AbstractVector, zs::AbstractMatrix) -> nodes, faces\n\nConvert a regular grid (x, y) and corresponding z-coordinates to triangular facets\n\n| ⧹| ⧹| ⧹|\n\nj+1 ・–C–D–・     |⧹ |⧹ |⧹ |     | ⧹| ⧹| ⧹| j   ・–A–B–・     |⧹ |⧹ |⧹ |        i  i+1\n\nArguments\n\nxs::AbstractVector : x-coordinates of grid points (should be sorted)\nys::AbstractVector : y-coordinates of grid points (should be sorted)\nzs::AbstractMatrix : z-coordinates of grid points\n\n\n\n\n\n","category":"method"},{"location":"api/#AsteroidThermoPhysicalModels.init_temperature!-Tuple{AsteroidThermoPhysicalModels.BinaryTPM, Real}","page":"API","title":"AsteroidThermoPhysicalModels.init_temperature!","text":"init_temperature!(btpm::BinaryTPM, T₀::Real)\n\nInitialize all temperature cells at the given temperature T₀\n\nArguments\n\nbtpm : Thermophysical model for a binary asteroid\nT₀   : Initial temperature of all cells [K]\n\n\n\n\n\n","category":"method"},{"location":"api/#AsteroidThermoPhysicalModels.init_temperature!-Tuple{AsteroidThermoPhysicalModels.SingleTPM, Real}","page":"API","title":"AsteroidThermoPhysicalModels.init_temperature!","text":"init_temperature!(stpm::SingleTPM, T₀::Real)\n\nInitialize all temperature cells at the given temperature T₀\n\nArguments\n\nstpm : Thermophysical model for a single asteroid\nT₀   : Initial temperature of all cells [K]\n\n\n\n\n\n","category":"method"},{"location":"api/#AsteroidThermoPhysicalModels.isilluminated-Tuple{ShapeModel, StaticArraysCore.StaticArray{Tuple{3}, T, 1} where T, Integer}","page":"API","title":"AsteroidThermoPhysicalModels.isilluminated","text":"isilluminated(shape::ShapeModel, r☉::StaticVector{3}, i::Integer) -> Bool\n\nReturn if the i-th face of the shape model is illuminated by the direct sunlight or not\n\nArguments\n\nshape : Shape model of an asteroid\nr☉    : Sun's position in the asteroid-fixed frame, which doesn't have to be normalized.\ni     : Index of the face to be checked\n\n\n\n\n\n","category":"method"},{"location":"api/#AsteroidThermoPhysicalModels.load_shape_grid-Tuple{AbstractVector, AbstractVector, AbstractMatrix}","page":"API","title":"AsteroidThermoPhysicalModels.load_shape_grid","text":"load_shape_grid(xs, ys, zs; scale=1.0, find_visible_facets=false) -> shape\n\nConvert a regular grid (x, y) to a shape model\n\nArguments\n\nxs::AbstractVector : x-coordinates of grid points\nys::AbstractVector : y-coordinates of grid points\nzs::AbstractMatrix : z-coordinates of grid points\n\n\n\n\n\n","category":"method"},{"location":"api/#AsteroidThermoPhysicalModels.loadobj-Tuple{String}","page":"API","title":"AsteroidThermoPhysicalModels.loadobj","text":"loadobj(shapepath::String; scale=1, message=true) -> nodes, faces\n\n\n\n\n\n","category":"method"},{"location":"api/#AsteroidThermoPhysicalModels.mutual_heating!-Tuple{AsteroidThermoPhysicalModels.BinaryTPM, Any, Any}","page":"API","title":"AsteroidThermoPhysicalModels.mutual_heating!","text":"mutual_heating!(btpm::BinaryTPM, rₛ, R₂₁)\n\nCalculate the mutual heating between the primary and secondary asteroids.\n\nArguments\n\nbtpm : Thermophysical model for a binary asteroid\nrₛ   : Position of the secondary relative to the primary (NOT normalized)\nR₂₁  : Rotation matrix from secondary to primary\n\nTO DO\n\nNeed to consider local horizon?\n\n\n\n\n\n","category":"method"},{"location":"api/#AsteroidThermoPhysicalModels.mutual_shadowing!-Tuple{AsteroidThermoPhysicalModels.BinaryTPM, Any, Any, Any}","page":"API","title":"AsteroidThermoPhysicalModels.mutual_shadowing!","text":"mutual_shadowing!(btpm::BinaryTPM, r☉, rₛ, R₂₁)\n\nDetect eclipse events between the primary and secondary, and update the solar fluxes of the faces.\n\nArguments\n\nbtpm : Thermophysical model for a binary asteroid\nr☉   : Position of the sun relative to the primary       (NOT normalized)\nrₛ   : Position of the secondary relative to the primary (NOT normalized)\nR₂₁  : Rotation matrix from secondary to primary\n\n\n\n\n\n","category":"method"},{"location":"api/#AsteroidThermoPhysicalModels.polyhedron_volume-Tuple{Any, Any}","page":"API","title":"AsteroidThermoPhysicalModels.polyhedron_volume","text":"polyhedron_volume(nodes, faces)      -> vol\npolyhedron_volume(shape::ShapeModel) -> vol\n\nCalculate volume of a polyhedral\n\n\n\n\n\n","category":"method"},{"location":"api/#AsteroidThermoPhysicalModels.raycast-NTuple{4, StaticArraysCore.StaticArray{Tuple{3}, T, 1} where T}","page":"API","title":"AsteroidThermoPhysicalModels.raycast","text":"raycast(A, B, C, R) -> Bool\n\nIntersection detection between ray R and triangle ABC. Note that the starting point of the ray is the origin (0, 0, 0).\n\n\n\n\n\n","category":"method"},{"location":"api/#AsteroidThermoPhysicalModels.raycast-NTuple{5, StaticArraysCore.StaticArray{Tuple{3}, T, 1} where T}","page":"API","title":"AsteroidThermoPhysicalModels.raycast","text":"raycast(A, B, C, R, O) -> Bool\n\nIntersection detection between ray R and triangle ABC. Use when the starting point of the ray is an arbitrary point O.\n\n\n\n\n\n","category":"method"},{"location":"api/#AsteroidThermoPhysicalModels.run_TPM!-Tuple{AsteroidThermoPhysicalModels.BinaryTPM, Any, Vector{Float64}, Vector{Int64}, Vector{Int64}}","page":"API","title":"AsteroidThermoPhysicalModels.run_TPM!","text":"run_TPM!(btpm::BinaryTPM, ephem, savepath)\n\nRun TPM for a binary asteroid.\n\nArguments\n\nbtpm          : Thermophysical model for a binary asteroid\nephem         : Ephemerides\ntime : Ephemeris times\nsun1 : Sun's position in the primary's frame\nsun2 : Sun's position in the secondary's frame\nsec  : Secondary's position in the primary's frame\nP2S  : Rotation matrix from primary to secondary frames\nS2P  : Rotation matrix from secondary to primary frames\ntimes_to_save : Timesteps to save temperature\nface_ID_pri   : Face indices where to save subsurface termperature for the primary\nface_ID_sec   : Face indices where to save subsurface termperature for the secondary\n\nKeyword arguments\n\nshow_progress : Flag to show the progress meter\n\n\n\n\n\n","category":"method"},{"location":"api/#AsteroidThermoPhysicalModels.run_TPM!-Tuple{AsteroidThermoPhysicalModels.SingleTPM, Any, Vector{Float64}, Vector{Int64}}","page":"API","title":"AsteroidThermoPhysicalModels.run_TPM!","text":"run_TPM!(stpm::SingleTPM, ephem, savepath)\n\nRun TPM for a single asteroid.\n\nArguments\n\nstpm          : Thermophysical model for a single asteroid\nephem         : Ephemerides\nephem.time : Ephemeris times\nephem.sun  : Sun's position in the asteroid-fixed frame (Not normalized)\ntimes_to_save : Timesteps to save temperature\nface_ID       : Face indices where to save subsurface termperature\n\nKeyword arguments\n\nshow_progress : Flag to show the progress meter\n\n\n\n\n\n","category":"method"},{"location":"api/#AsteroidThermoPhysicalModels.subsolar_temperature-Tuple{Any, AsteroidThermoPhysicalModels.AbstractThermoParams}","page":"API","title":"AsteroidThermoPhysicalModels.subsolar_temperature","text":"subsolar_temperature(r☉) -> Tₛₛ\n\nSubsolar temperature [K] on an asteroid at a heliocentric distance r☉ [m], assuming radiative equilibrium with zero conductivity.\n\n\n\n\n\n","category":"method"},{"location":"api/#AsteroidThermoPhysicalModels.surface_temperature-Tuple{AsteroidThermoPhysicalModels.SingleTPM}","page":"API","title":"AsteroidThermoPhysicalModels.surface_temperature","text":"Return surface temperature of a single asteroid corrsponding to each face.\n\n\n\n\n\n","category":"method"},{"location":"api/#AsteroidThermoPhysicalModels.thermal_inertia-Tuple{Any, Any, Any}","page":"API","title":"AsteroidThermoPhysicalModels.thermal_inertia","text":"thermal_inertia(k, ρ, Cp) -> Γ\n\nArguments\n\nk  : Thermal conductivity [W/m/K]\nρ  : Material density [kg/m³]\nCₚ : Heat capacity [J/kg/K]\n\nReturn\n\nΓ : Thermal inertia [J ⋅ m⁻² ⋅ K⁻¹ ⋅ s⁻⁰⁵ (tiu)]\n\n\n\n\n\n","category":"method"},{"location":"api/#AsteroidThermoPhysicalModels.thermal_skin_depth-NTuple{4, Any}","page":"API","title":"AsteroidThermoPhysicalModels.thermal_skin_depth","text":"thermal_skin_depth(P, k, ρ, Cp) -> l_2π\n\nArguments\n\nP  : Cycle of thermal cycle [sec]\nk  : Thermal conductivity [W/m/K]\nρ  : Material density [kg/m³]\nCₚ : Heat capacity [J/kg/K]\n\nReturn\n\nl_2π : Thermal skin depth [m], as defined in Rozitis & Green (2011).\n\n\n\n\n\n","category":"method"},{"location":"api/#AsteroidThermoPhysicalModels.thermoparams-Tuple{}","page":"API","title":"AsteroidThermoPhysicalModels.thermoparams","text":"thermoparams(; A_B, A_TH, k, ρ, Cp, ε, t_begin, t_end, Nt, z_max, Nz, P)\n\n\n\n\n\n","category":"method"},{"location":"api/#AsteroidThermoPhysicalModels.tridiagonal_matrix_algorithm!-NTuple{5, Any}","page":"API","title":"AsteroidThermoPhysicalModels.tridiagonal_matrix_algorithm!","text":"tridiagonal_matrix_algorithm!(a, b, c, d, x)\ntridiagonal_matrix_algorithm!(stpm::SingleTPM)\n\nTridiagonal matrix algorithm to solve the heat conduction equation by the backward Euler and Crank-Nicolson methods.\n\n| b₁ c₁ 0  ⋯  0   | | x₁ |   | d₁ |\n| a₂ b₂ c₂ ⋯  0   | | x₂ |   | d₂ |\n| 0  a₃ b₃ ⋯  0   | | x₃ | = | d₃ |\n| ⋮  ⋮  ⋮  ⋱  cₙ₋₁| | ⋮  |   | ⋮  |\n| 0  0  0  aₙ bₙ  | | xₙ |   | dₙ |\n\nReferences\n\nhttps://en.wikipedia.org/wiki/Tridiagonalmatrixalgorithm\n\n\n\n\n\n","category":"method"},{"location":"api/#AsteroidThermoPhysicalModels.update_TPM_result!-Tuple{AsteroidThermoPhysicalModels.BinaryTPMResult, AsteroidThermoPhysicalModels.BinaryTPM, Integer}","page":"API","title":"AsteroidThermoPhysicalModels.update_TPM_result!","text":"update_TPM_result!(result::BinaryTPMResult, btpm::BinaryTPM, ephem, nₜ::Integer)\n\nSave the results of TPM at the time step nₜ to result.\n\nArguments\n\nresult : Output data format for BinaryTPM\nbtpm   : Thermophysical model for a binary asteroid\nephem  : Ephemerides\nnₜ     : Time step\n\n\n\n\n\n","category":"method"},{"location":"api/#AsteroidThermoPhysicalModels.update_TPM_result!-Tuple{AsteroidThermoPhysicalModels.SingleTPMResult, AsteroidThermoPhysicalModels.SingleTPM, Integer}","page":"API","title":"AsteroidThermoPhysicalModels.update_TPM_result!","text":"update_TPM_result!(result::SingleTPMResult, stpm::SingleTPM, nₜ::Integer)\n\nSave the results of TPM at the time step nₜ to result.\n\nArguments\n\nresult : Output data format for SingleTPM\nstpm   : Thermophysical model for a single asteroid\nnₜ     : Time step to save data\n\n\n\n\n\n","category":"method"},{"location":"api/#AsteroidThermoPhysicalModels.update_flux_rad_single!-Tuple{AsteroidThermoPhysicalModels.BinaryTPM}","page":"API","title":"AsteroidThermoPhysicalModels.update_flux_rad_single!","text":"update_flux_rad_single!(btpm::BinaryTPM)\n\nUpdate flux of absorption of thermal radiation from surrounding surface. Single radiation-absorption is only considered, assuming albedo is close to zero at thermal infrared wavelength.\n\nArguments\n\nbtpm : Thermophysical model for a binary asteroid\n\n\n\n\n\n","category":"method"},{"location":"api/#AsteroidThermoPhysicalModels.update_flux_rad_single!-Tuple{AsteroidThermoPhysicalModels.SingleTPM}","page":"API","title":"AsteroidThermoPhysicalModels.update_flux_rad_single!","text":"update_flux_rad_single!(stpm::SingleTPM)\n\nUpdate flux of absorption of thermal radiation from surrounding surface. Single radiation-absorption is only considered, assuming albedo is close to zero at thermal infrared wavelength.\n\nArguments\n\nstpm : Thermophysical model for a single asteroid\n\n\n\n\n\n","category":"method"},{"location":"api/#AsteroidThermoPhysicalModels.update_flux_scat_single!-Tuple{AsteroidThermoPhysicalModels.BinaryTPM}","page":"API","title":"AsteroidThermoPhysicalModels.update_flux_scat_single!","text":"update_flux_scat_single!(btpm::BinaryTPM)\n\nUpdate flux of scattered sunlight, only considering single scattering.\n\nArguments\n\nbtpm : Thermophysical model for a binary asteroid\n\n\n\n\n\n","category":"method"},{"location":"api/#AsteroidThermoPhysicalModels.update_flux_scat_single!-Tuple{AsteroidThermoPhysicalModels.SingleTPM}","page":"API","title":"AsteroidThermoPhysicalModels.update_flux_scat_single!","text":"update_flux_scat_single!(stpm::SingleTPM)\n\nUpdate flux of scattered sunlight, only considering single scattering.\n\nArguments\n\nstpm : Thermophysical model for a single asteroid\n\n\n\n\n\n","category":"method"},{"location":"api/#AsteroidThermoPhysicalModels.update_flux_sun!-Tuple{AsteroidThermoPhysicalModels.BinaryTPM, StaticArraysCore.StaticArray{Tuple{3}, T, 1} where T, StaticArraysCore.StaticArray{Tuple{3}, T, 1} where T}","page":"API","title":"AsteroidThermoPhysicalModels.update_flux_sun!","text":"update_flux_sun!(btpm::BinaryTPM, r☉₁::StaticVector{3}, r☉₂::StaticVector{3})\n\nArguments\n\nbtpm : Thermophysical model for a binary asteroid\nr☉₁  : Sun's position in the body-fixed frame of the primary, which is not normalized.\nr☉₂  : Sun's position in the body-fixed frame of the secondary, which is not normalized.\n\n\n\n\n\n","category":"method"},{"location":"api/#AsteroidThermoPhysicalModels.update_flux_sun!-Tuple{AsteroidThermoPhysicalModels.SingleTPM, StaticArraysCore.StaticArray{Tuple{3}, T, 1} where T, Real}","page":"API","title":"AsteroidThermoPhysicalModels.update_flux_sun!","text":"update_flux_sun!(stpm::SingleTPM, r̂☉::StaticVector{3}, F☉::Real)\n\nUpdate solar irradiation flux on every face of a shape model.\n\nshape : Shape model\nr̂☉    : Normalized vector indicating the direction of the sun in the body-fixed frame\nF☉    : Solar radiation flux [W/m²]\n\n\n\n\n\n","category":"method"},{"location":"api/#AsteroidThermoPhysicalModels.update_flux_sun!-Tuple{AsteroidThermoPhysicalModels.SingleTPM, StaticArraysCore.StaticArray{Tuple{3}, T, 1} where T}","page":"API","title":"AsteroidThermoPhysicalModels.update_flux_sun!","text":"update_flux_sun!(stpm::SingleTPM, r☉::StaticVector{3})\n\nUpdate solar irradiation flux on every face of a shape model.\n\nArguments\n\nstpm : Thermophysical model for a single asteroid\nr☉   : Position of the sun in the body-fixed frame (NOT normalized)\n\n\n\n\n\n","category":"method"},{"location":"api/#AsteroidThermoPhysicalModels.update_lower_temperature!-Tuple{AsteroidThermoPhysicalModels.SingleTPM}","page":"API","title":"AsteroidThermoPhysicalModels.update_lower_temperature!","text":"update_bottom_temperature!(shape::ShapeModel)\n\nUpdate the temperature of the bottom surface based on the boundary condition stpm.BC_LOWER.\n\nArguments\n\nstpm       : Thermophysical model for a single asteroid\n\n\n\n\n\n","category":"method"},{"location":"api/#AsteroidThermoPhysicalModels.update_surface_temperature!-Tuple{AbstractVector, Vararg{Float64, 6}}","page":"API","title":"AsteroidThermoPhysicalModels.update_surface_temperature!","text":"update_surface_temperature!(T::AbstractVector, F_total::Real, k::Real, l::Real, Δz::Real, ε::Real)\n\nNewton's method to update the surface temperature under radiation boundary condition.\n\nArguments\n\nT       : 1-D array of temperatures\nF_total : Total energy absorbed by the facet\nΓ       : Thermal inertia [tiu]\nP       : Period of thermal cycle [sec]\nΔz̄      : Non-dimensional step in depth, normalized by thermal skin depth l\nε       : Emissivity\n\n\n\n\n\n","category":"method"},{"location":"api/#AsteroidThermoPhysicalModels.update_temperature!-Tuple{AsteroidThermoPhysicalModels.BinaryTPM, Any}","page":"API","title":"AsteroidThermoPhysicalModels.update_temperature!","text":"update_temperature!(btpm::BinaryTPM, Δt)\n\nCalculate the temperature for the next time step based on 1D heat conductivity equation.\n\nArguments\n\nbtpm : Thermophysical model for a binary asteroid\nΔt   : Time step [sec]\n\n\n\n\n\n","category":"method"},{"location":"api/#AsteroidThermoPhysicalModels.update_temperature!-Tuple{AsteroidThermoPhysicalModels.SingleTPM, Any}","page":"API","title":"AsteroidThermoPhysicalModels.update_temperature!","text":"update_temperature!(stpm::SingleTPM, Δt)\n\nCalculate the temperature for the next time step based on 1D heat conduction equation.\n\nArguments\n\nstpm : Thermophysical model for a single asteroid\nΔt   : Time step [sec]\n\n\n\n\n\n","category":"method"},{"location":"api/#AsteroidThermoPhysicalModels.update_thermal_force!-Tuple{AsteroidThermoPhysicalModels.BinaryTPM}","page":"API","title":"AsteroidThermoPhysicalModels.update_thermal_force!","text":"update_thermal_force!(btpm::BinaryTPM)\n\nCalculate the thermal force and torque on every face and integrate them over all faces.\n\nArguments\n\nbtpm : Thermophysical model for a binary asteroid\n\n\n\n\n\n","category":"method"},{"location":"api/#AsteroidThermoPhysicalModels.update_thermal_force!-Tuple{AsteroidThermoPhysicalModels.SingleTPM}","page":"API","title":"AsteroidThermoPhysicalModels.update_thermal_force!","text":"update_thermal_force!(stpm::SingleTPM)\n\nCalculate the thermal force and torque on every face and integrate them over all faces.\n\nArguments\n\nstpm : Thermophysical model for a single asteroid\n\n\n\n\n\n","category":"method"},{"location":"api/#AsteroidThermoPhysicalModels.update_upper_temperature!-Tuple{AsteroidThermoPhysicalModels.SingleTPM, Integer}","page":"API","title":"AsteroidThermoPhysicalModels.update_upper_temperature!","text":"update_upper_temperature!(stpm::SingleTPM, nₛ::Integer)\n\nUpdate the temperature of the upper surface based on the boundary condition stpm.BC_UPPER.\n\nArguments\n\nstpm      : Thermophysical model for a single asteroid\nnₛ        : Index of the face of the shape model\n\n\n\n\n\n","category":"method"},{"location":"api/#AsteroidThermoPhysicalModels.view_factor-NTuple{5, Any}","page":"API","title":"AsteroidThermoPhysicalModels.view_factor","text":"view_factor(cᵢ, cⱼ, n̂ᵢ, n̂ⱼ, aⱼ) -> fᵢⱼ, dᵢⱼ, d̂ᵢⱼ\n\nView factor from facet i to j, assuming Lambertian emission.\n\n\n(i)   fᵢⱼ   (j)\n△    –>    △\n\ncᵢ          cⱼ  : Center of each face\nn̂ᵢ          n̂ⱼ  : Normal vector of each face\n      aⱼ  : Area of j-th face\n\n\n\n\n\n\n","category":"method"},{"location":"physical_model/#Physical-model","page":"Physical model","title":"Physical model","text":"","category":"section"},{"location":"physical_model/#Overview","page":"Physical model","title":"Overview","text":"","category":"section"},{"location":"physical_model/#Symbols","page":"Physical model","title":"Symbols","text":"","category":"section"},{"location":"physical_model/","page":"Physical model","title":"Physical model","text":"Symbol Unit Description\nt [s] time\nT [K] Termperture\nA_B [-] Bond albedo\nA_textTH [-] Albedo at thermal infrared wavelength\nF_textsun [W/m²] Flux of direct sunlight\nF_textscat [W/m²] Flux of scattered light\nF_textrad [W/m²] Flux of thermal radiation from surrounding surface\nrho [kg/m³] Density\nC_p [J/K] Heat capacity at constant pressure\nP [s] Rotation period\nl [m] Thermal skin depth\nk [W/m/K] Thermal conductivity\nz [m] Depth\nE [J] Emittance energy\nGamma [J ⋅ m⁻² ⋅ K⁻¹ ⋅ s⁻⁰⁵ (tiu)] Thermal inertia\nvarepsilon [-] Emissivity\nPhi [W/m²] Solar energy flux","category":"page"},{"location":"#AsteroidThermoPhysicalModels.jl","page":"Home","title":"AsteroidThermoPhysicalModels.jl","text":"","category":"section"},{"location":"","page":"Home","title":"Home","text":"A package for dynamical simulation of an asteroid.","category":"page"}]
}
