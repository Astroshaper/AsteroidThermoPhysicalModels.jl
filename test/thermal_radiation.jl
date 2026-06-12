#=
thermal_radiation.jl

Tests for thermal radiation calculations.
This test validates:
- Blackbody radiation calculations (Planck's law)
- Stefan-Boltzmann law implementation
- View factor calculations for thermal radiation exchange
- Self-heating effects between asteroid surface elements
=#

@testset "thermal_radiation" begin

    msg = """\n
    ┏━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━┓
    |                Test: thermal_radiation                 |
    ┗━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━┛
    """
    println(msg)

    ## Calculate the intensity of blackbody radiation of at a wavelength of 6e-7 m and a temperature of 5850 K,
    ## and compare it with the value calculated by Planck.jl.
    ## cf. https://github.com/JuliaAstro/Planck.jl/blob/main/src/Planck.jl
    @test AsteroidThermoPhysicalModels.blackbody_radiance(6e-7, 5850) ≈ 2.583616647617974e13

    ## Check the value of Stefan-Boltzmann law at a temperature of 5850 K
    @test AsteroidThermoPhysicalModels.blackbody_radiance(5850) ≈ AsteroidThermoPhysicalModels.σ_SB * 5850^4

    ######## Thermal radiance from the local terrain model ########

    ## --- Load shape model ---
    path_obj = joinpath("shape", "fractal_v2572_f5000.obj")
    shape = load_shape_obj(path_obj; scale=1, with_face_visibility=true)
    n_face = length(shape.faces)  # Number of faces

    ## --- Ephemerides ---
    P = SPICE.convrt(8, "hours", "seconds")  # Rotation period of the asteroid [s]
    
    ncycles = 2  # Number of cycles to perform TPM
    nsteps_in_cycle = 360  # Number of time steps in one rotation period
    
    ## TPM simulation time duration (ephemerides time)
    et_begin = 0.0  # Start time of TPM
    et_end   = et_begin + P * ncycles  # End time of TPM
    et_range = range(et_begin, et_end; length=nsteps_in_cycle*ncycles+1)

    ## Rotation matrix between the global (asteroid-fixed) frame and the local frame at `(lat, lon)`.
    local_to_global(lat, lon) = RotZYZ(lon, π/2 - lat, π/2)
    global_to_local(lat, lon) = inv(local_to_global(lat, lon))
    
    ## Solar position in the local frame of the roughness map
    lat = deg2rad(0)  # Latitude of the local frame [rad]
    lon = deg2rad(0)  # Longitude of the local frame [rad]

    r☉₀ = [SPICE.convrt(-1, "au", "m"), 0, 0]  # Sun's position at the initial time step in the asteoid-fixed frame

    times = collect(et_range)
    r_sun = [global_to_local(lat, lon) * inv(RotZ(2π * et / P)) * r☉₀ for et in et_range]  # Sun's position in local frame [m]

    ephem = SingleAsteroidEphemerides(times, r_sun)
    
    ## --- Thermal properties ---
    k  = 0.1     # Thermal conductivity [W/m/K]
    ρ  = 1270.0  # Density [kg/m³]
    Cₚ = 600.0   # Heat capacity [J/kg/K]
    
    l = AsteroidThermoPhysicalModels.thermal_skin_depth(P, k, ρ, Cₚ)
    Γ = AsteroidThermoPhysicalModels.thermal_inertia(k, ρ, Cₚ)

    R_vis = 0.04  # Reflectance in visible light [-]
    R_ir  = 0.0   # Reflectance in thermal infrared [-]
    ε     = 1.0   # Emissivity [-]

    z_max = 0.6   # Depth of the lower boundary of a heat conduction equation [m]
    n_depth = 41  # Number of depth steps
    Δz = z_max / (n_depth - 1)  # Depth step width [m]

    thermo_params = ThermoParams(k, ρ, Cₚ, R_vis, R_ir, ε, z_max, Δz, n_depth)

    ## --- Setting of TPM ---
    problem = SingleAsteroidThermoPhysicalProblem(shape, thermo_params;
        with_self_shadowing      = true,
        with_self_heating        = true,
        upper_boundary_condition = RadiationBoundaryCondition(),
        lower_boundary_condition = InsulationBoundaryCondition(),
    )

    times_to_save = ephem.times[end-nsteps_in_cycle:end]  # Save temperature during the final rotation
    face_ID = [49, 340, 648]  # Face indices to save subsurface temperature

    solution = solve(problem, ExplicitEuler();
        ephem         = ephem,
        times_to_save = times_to_save,
        face_ID       = face_ID,
        T₀            = 200.0,
    )

    ## --- Check the thermal radiation from the local terrain model ---
    obs_above = SVector{3, Float64}(0, 0, 1000)  # Observer is just above the local terrain model
    obs_east  = RotY(+π/6) * obs_above           # Observed from 30° east
    obs_west  = RotY(-π/6) * obs_above           # Observed from 30° west

    emissivities = fill(1.0, length(shape.faces))
    temperatures = solution.surface_temperature[:, 181]

    ## Expected values of the thermal radiation
    ## - Observation from 30° east   : 19.89394442347112  [W/m²]
    ## - Observation from just above : 18.01010231251351  [W/m²]
    ## - Observation from 30° west   : 12.406927167050457 [W/m²]
    @test AsteroidThermoPhysicalModels.thermal_radiance(shape, emissivities, temperatures, obs_east)  ≈ 19.89394442347112
    @test AsteroidThermoPhysicalModels.thermal_radiance(shape, emissivities, temperatures, obs_above) ≈ 18.01010231251351
    @test AsteroidThermoPhysicalModels.thermal_radiance(shape, emissivities, temperatures, obs_west)  ≈ 12.406927167050457
end
