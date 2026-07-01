#=
TPM_zero-conductivity.jl

Tests the special case of zero thermal conductivity.
This test validates:
- Instantaneous thermal equilibrium (no heat conduction)
- Surface temperatures determined purely by energy balance
- Limiting case behavior for very low thermal inertia materials
=#

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

    ## --- Ephemerides ---
    P = SPICE.convrt(8, "hours", "seconds")  # Rotation period of the asteroid [s]
    n_cycle = 2                              # Number of cycles to perform TPM
    n_step_in_cycle = 72                     # Number of time steps in one rotation period

    et_begin = 0.0                     # Start time of TPM
    et_end   = et_begin + P * n_cycle  # End time of TPM
    et_range = range(et_begin, et_end; length=n_step_in_cycle*n_cycle+1)

    times = collect(et_range)
    r_sun = [inv(RotZ(2π * et / P)) * [SPICE.convrt(1, "au", "m"), 0, 0] for et in et_range]  # Sun's position in asteroid-fixed frame [m]

    ephem = SingleAsteroidEphemerides(times, r_sun)

    ## --- Load obj file ---
    path_obj = joinpath("shape", "ryugu_test.obj")
    shape = load_shape_obj(path_obj; scale=1000, with_face_visibility=true, with_bvh=true)
    n_face = length(shape.faces)  # Number of faces

    ## --- Thermal properties: zero-conductivity case ---
    thermo_params = ThermoParams(
        conductivity    = 0.0,    # Thermal conductivity [W/m/K]
        density         = 1270.0, # Density [kg/m³]
        heat_capacity   = 600.0,  # Heat capacity [J/kg/K]
        reflectance_vis = 0.1,    # Reflectance in visible light [-]
        reflectance_ir  = 0.0,    # Reflectance in thermal infrared [-]
        emissivity      = 1.0,    # Emissivity [-]
    )
    grid_params = GridParams(; z_max=0.6, n_depth=41)

    ## --- Setting of TPM ---
    problem = SingleAsteroidThermoPhysicalProblem(shape, thermo_params, grid_params;
        with_self_shadowing      = true,
        with_self_heating        = false,
        upper_boundary_condition = RadiationBoundaryCondition(),
        lower_boundary_condition = InsulationBoundaryCondition(),
    )

    ## --- Run TPM ---
    output_times        = ephem.times[end-n_step_in_cycle:end]  # Save temperature during the final rotation
    subsurface_face_ids = [1, 2]  # Face indices to save subsurface temperature
    output = SingleAsteroidOutputSpec(output_times, subsurface_face_ids)

    solution = solve(problem, ExplicitEuler();
        ephem               = ephem,
        output              = output,
        initial_temperature = 0.0,
    )

    ## --- Save TPM result ---
    @testset "Save TPM result" begin
        export_solution(DIR_OUTPUT, solution)

        @test isfile(joinpath(DIR_OUTPUT, "diagnostics.csv"))
        @test isfile(joinpath(DIR_OUTPUT, "surface_temperature.csv"))
        @test isfile(joinpath(DIR_OUTPUT, "subsurface_temperature.csv"))
        @test !isfile(joinpath(DIR_OUTPUT, "thermal_face_forces.csv"))  # save_face_forces=false by default
    end
end
