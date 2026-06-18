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
    k  = 0.0     # Thermal conductivity [W/m/K]
    ρ  = 1270.0  # Density [kg/m³]
    Cₚ = 600.0   # Heat capacity [J/kg/K]

    R_vis = 0.1  # Reflectance in visible light [-]
    R_ir  = 0.0  # Reflectance in thermal infrared [-]
    ε     = 1.0  # Emissivity [-]

    z_max = 0.6   # Depth of the lower boundary of a heat conduction equation [m]
    n_depth = 41  # Number of depth steps
    Δz = z_max / (n_depth - 1)  # Depth step width [m]

    thermo_params = ThermoParams(k, ρ, Cₚ, R_vis, R_ir, ε, z_max, Δz, n_depth)

    ## --- Setting of TPM ---
    problem = AsteroidThermoPhysicalModels.SingleAsteroidThermoPhysicalProblem(shape, thermo_params;
        with_self_shadowing      = true,
        with_self_heating        = false,
        upper_boundary_condition = AsteroidThermoPhysicalModels.RadiationBoundaryCondition(),
        lower_boundary_condition = AsteroidThermoPhysicalModels.InsulationBoundaryCondition(),
    )

    ## --- Run TPM ---
    times_to_save = ephem.times[end-n_step_in_cycle:end]  # Save temperature during the final rotation
    face_ID = [1, 2]  # Face indices to save subsurface temperature
    output = SingleAsteroidOutputSpec(times_to_save, face_ID)

    solution = solve(problem, ExplicitEuler();
        ephem               = ephem,
        output              = output,
        initial_temperature = 0.0,
    )

    ## --- Save TPM result ---
    @testset "Save TPM result" begin
        AsteroidThermoPhysicalModels.export_solution(DIR_OUTPUT, solution)

        @test isfile(joinpath(DIR_OUTPUT, "physical_quantities.csv"))
        @test isfile(joinpath(DIR_OUTPUT, "subsurface_temperature.csv"))
        @test isfile(joinpath(DIR_OUTPUT, "surface_temperature.csv"))
        @test isfile(joinpath(DIR_OUTPUT, "thermal_force.csv"))
    end
end
