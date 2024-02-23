# See https://github.com/Astroshaper/Astroshaper-examples/tree/main/TPM_Ryugu for more information.
@testset "TPM_zero-conductivity" begin
    DIR_OUTPUT = joinpath(@__DIR__, "output")
    rm(DIR_OUTPUT; recursive=true, force=true)
    mkpath(DIR_OUTPUT)

    msg = """\n
    ┏━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━┓
    |              Test: TPM_zero-conductivity               |
    ┗━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━┛
    """
    println(msg)

    ##= Ephemerides =##
    P = SPICE.convrt(8, "hours", "seconds")  # Rotation period of the asteroid [s]
    ncycles = 2                              # Number of cycles to perform TPM
    nsteps_in_cycle = 72                     # Number of time steps in one rotation period

    et_begin = 0.0                     # Start time of TPM
    et_end   = et_begin + P * ncycles  # End time of TPM
    et_range = range(et_begin, et_end; length=nsteps_in_cycle*ncycles+1)

    """
    - `time` : Ephemeris times
    - `sun`  : Sun's position in the asteroid-fixed frame
    """
    ephem = (
        time = collect(et_range),
        sun  = [inv(RotZ(2π * et / P)) * [SPICE.convrt(1, "au", "m"), 0, 0] for et in et_range],
    )

    ##= Load obj file =##
    path_obj = joinpath("shape", "ryugu_test.obj")
    shape = AsteroidThermoPhysicalModels.load_shape_obj(path_obj; scale=1000, find_visible_facets=true)

    ##= Thermal properties: zero-conductivity case =##
    thermo_params = AsteroidThermoPhysicalModels.thermoparams(
        P       = P,
        l       = 0.0,
        Γ       = 0.0,
        A_B     = 0.1,
        A_TH    = 0.0,
        ε       = 1.0,
        z_max   = 0.6,
        n_depth = 11,
    )

    println(thermo_params)

    ##= Setting of TPM =##
    stpm = AsteroidThermoPhysicalModels.SingleTPM(shape, thermo_params;
        SELF_SHADOWING = true,
        SELF_HEATING   = false,
        SOLVER         = AsteroidThermoPhysicalModels.ForwardEulerSolver(thermo_params),
        BC_UPPER       = AsteroidThermoPhysicalModels.RadiationBoundaryCondition(),
        BC_LOWER       = AsteroidThermoPhysicalModels.InsulationBoundaryCondition(),
    )
    AsteroidThermoPhysicalModels.init_temperature!(stpm, 0)

    ##= Run TPM =##
    times_to_save = ephem.time[end-nsteps_in_cycle:end]  # Save temperature during the final rotation
    face_ID = [1, 2]  # Face indices to save subsurface temperature

    result = AsteroidThermoPhysicalModels.run_TPM!(stpm, ephem, times_to_save, face_ID)
    
    ##= Save TPM result =##
    @testset "Save TPM result" begin
        AsteroidThermoPhysicalModels.export_TPM_results(DIR_OUTPUT, result)

        @test isfile(joinpath(DIR_OUTPUT, "physical_quantities.csv"))
        @test isfile(joinpath(DIR_OUTPUT, "subsurface_temperature.csv"))
        @test isfile(joinpath(DIR_OUTPUT, "surface_temperature.csv"))
        @test isfile(joinpath(DIR_OUTPUT, "thermal_force.csv"))
    end
end
