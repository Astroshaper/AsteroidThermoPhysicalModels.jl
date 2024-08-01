@testset "blackbody_radiation" begin

    msg = """\n
    ┏━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━┓
    |               Test: blackbody_radiation                |
    ┗━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━┛
    """
    println(msg)

    @test AsteroidThermoPhysicalModels.blackbody_radiation(6e-7, 5850) ≈ 2.583616647617974e13

    ##= Load shape model =##
    path_obj = joinpath("shape", "fractal_v2572_f5000.obj")
    shape = AsteroidThermoPhysicalModels.load_shape_obj(path_obj; scale=1, find_visible_facets=true)

    ##= Ephemerides =##
    P = SPICE.convrt(8, "hours", "seconds")  # Rotation period of the asteroid [s]
    
    ncycles = 2  # Number of cycles to perform TPM
    nsteps_in_cycle = 360  # Number of time steps in one rotation period
    
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

    """
    - `time` : Ephemeris times
    - `sun`  : Sun's position in the local frame
    """
    ephem = (
        time = collect(et_range),
        sun  = [global_to_local(lat, lon) * inv(RotZ(2π * et / P)) * r☉₀ for et in et_range],
    )
    
    ##= Thermal properties =##
    k  = 0.1
    ρ  = 1270.0
    Cₚ = 600.0
    
    l = AsteroidThermoPhysicalModels.thermal_skin_depth(P, k, ρ, Cₚ)
    Γ = AsteroidThermoPhysicalModels.thermal_inertia(k, ρ, Cₚ)

    thermo_params = AsteroidThermoPhysicalModels.thermoparams(
        P       = P,
        l       = l,
        Γ       = Γ,
        R_vis   = 0.04,
        R_ir    = 0.0,
        ε       = 1.0,
        z_max   = 0.6,
        n_depth = 41,
    )

    println(thermo_params)

    ##= Setting of TPM =##
    stpm = AsteroidThermoPhysicalModels.SingleTPM(shape, thermo_params;
        SELF_SHADOWING = true,
        SELF_HEATING   = true,
        SOLVER         = AsteroidThermoPhysicalModels.ForwardEulerSolver(thermo_params),
        BC_UPPER       = AsteroidThermoPhysicalModels.RadiationBoundaryCondition(),
        BC_LOWER       = AsteroidThermoPhysicalModels.InsulationBoundaryCondition(),
    )
    AsteroidThermoPhysicalModels.init_temperature!(stpm, 200)

    times_to_save = ephem.time[end-nsteps_in_cycle:end]  # Save temperature during the final rotation
    face_ID = [49, 340, 648]  # Face indices to save subsurface temperature
    
    result = run_TPM!(stpm, ephem, times_to_save, face_ID)

    ##= Check the thermal radiation from the local map =##
    obs_above = [0.0, 0.0, 1000.0]      # Observer is just above the local map
    obs_east  = RotY(+π/6) * obs_above  # Observed from 30° east
    obs_west  = RotY(-π/6) * obs_above  # Observed from 30° west

    emissivities = fill(1.0, length(shape.faces))
    temperatures = result.surface_temperature[:, 181]

    println("==== Thermal radiation from the local map ====")
    println("- Observation from 30° east   : ", AsteroidThermoPhysicalModels.thermal_radiation(shape, emissivities, temperatures, obs_east))
    println("- Observation from just above : ", AsteroidThermoPhysicalModels.thermal_radiation(shape, emissivities, temperatures, obs_above))
    println("- Observation from 30° west   : ", AsteroidThermoPhysicalModels.thermal_radiation(shape, emissivities, temperatures, obs_west))
end