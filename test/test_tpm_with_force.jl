#=
test_tpm_with_force.jl

Tests for the force/torque-enabled simulation paths:
    SingleAsteroidEphemerides{<:AbstractVector} with save_forces/save_torques=true
    BinaryAsteroidEphemerides{<:AbstractVector} with save_forces/save_torques=true

Covers: _alloc_solution (with forces/torques), record_timestep! (with R),
        export_solution, _solve {<:AbstractVector} (single + binary).
=#

@testset "TPM with force/torque" begin
    msg = """
    ┏━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━┓
    |             Test: TPM with force/torque                |
    ┗━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━┛
    """
    println(msg)

    au2m = AsteroidThermoPhysicalModels.au2m

    # Small orbit: 2 rotations, 36 steps/cycle
    P            = SPICE.convrt(8, "hours", "seconds")
    n_cycle      = 2
    n_step_cycle = 36
    et_range     = range(0.0, P * n_cycle; length=n_step_cycle * n_cycle + 1)
    times        = collect(et_range)

    # Rotation matrices: body-fixed rotates with period P around z-axis
    R_b2i = [RotZ(2π * et / P) for et in et_range]

    # Sun position in body-fixed frame (inverse rotation of inertial position)
    r_sun = [inv(RotZ(2π * et / P)) * SVector(au2m, 0.0, 0.0) for et in et_range]

    path_obj      = joinpath("shape", "icosahedron.obj")
    thermo_params = ThermoParams(0.1, 1000.0, 600.0, 0.1, 0.0, 0.9, 0.1, 0.01, 5)
    output_times        = times[end-n_step_cycle:end]
    subsurface_face_ids = [1]

    @testset "Single asteroid — {<:AbstractVector}" begin
        DIR_OUTPUT = mktempdir()

        ephem = SingleAsteroidEphemerides(times, r_sun, R_b2i)
        @test ephem isa SingleAsteroidEphemerides{Vector{SMatrix{3,3,Float64,9}}}

        shape   = load_shape_obj(path_obj; scale=1000, with_face_visibility=false, with_bvh=false)
        problem = SingleAsteroidThermoPhysicalProblem(shape, thermo_params;
            with_self_shadowing = false,
            with_self_heating   = false,
        )

        output = SingleAsteroidOutputSpec(output_times, subsurface_face_ids, true, true, false, true, true)
        solution = solve(problem, CrankNicolson();
            ephem               = ephem,
            output              = output,
            initial_temperature = 200.0,
            show_progress       = true,
        )

        # Solution type
        @test solution isa SingleAsteroidThermoPhysicalSolution

        # forces/torques are populated at output_times only
        @test solution.forces  isa Vector{SVector{3,Float64}}
        @test solution.torques isa Vector{SVector{3,Float64}}
        @test length(solution.forces)  == length(output_times)
        @test length(solution.torques) == length(output_times)

        # export: thermal_net_forces.csv must contain force/torque columns
        export_solution(DIR_OUTPUT, solution)
        @test isfile(joinpath(DIR_OUTPUT, "thermal_net_forces.csv"))
        df = CSV.read(joinpath(DIR_OUTPUT, "thermal_net_forces.csv"), DataFrame)
        @test "force_x"  in names(df)
        @test "force_y"  in names(df)
        @test "force_z"  in names(df)
        @test "torque_x" in names(df)
        @test "torque_y" in names(df)
        @test "torque_z" in names(df)
    end

    @testset "Single asteroid — save_face_forces=true" begin
        DIR_OUTPUT = mktempdir()

        # face forces are in the body-fixed frame: no rotation matrices needed
        ephem_no_rot = SingleAsteroidEphemerides(times, r_sun)

        shape   = load_shape_obj(path_obj; scale=1000, with_face_visibility=false, with_bvh=false)
        problem = SingleAsteroidThermoPhysicalProblem(shape, thermo_params;
            with_self_shadowing = false,
            with_self_heating   = false,
        )
        n_face = length(shape.faces)

        output = SingleAsteroidOutputSpec(output_times, subsurface_face_ids;
            save_surface_temperature    = false,
            save_subsurface_temperature = false,
            save_face_forces            = true,
        )
        solution = solve(problem, CrankNicolson();
            ephem               = ephem_no_rot,
            output              = output,
            initial_temperature = 200.0,
            show_progress       = false,
        )

        @test solution.face_forces isa Matrix{SVector{3, Float64}}
        @test size(solution.face_forces) == (n_face, length(output_times))

        # export: thermal_face_forces.csv must be written with correct columns
        export_solution(DIR_OUTPUT, solution)
        @test isfile(joinpath(DIR_OUTPUT, "thermal_face_forces.csv"))
        df = CSV.read(joinpath(DIR_OUTPUT, "thermal_face_forces.csv"), DataFrame)
        @test "force_x" in names(df)
        @test "force_y" in names(df)
        @test "force_z" in names(df)
        @test nrow(df) == n_face * length(output_times)
    end

    @testset "Binary asteroid — {<:AbstractVector}" begin
        DIR_OUTPUT = mktempdir()

        r_secondary            = [[1e4, 0.0, 0.0]          for _ in et_range]
        R_primary_to_secondary = [Matrix{Float64}(I, 3, 3) for _ in et_range]
        R_primary_to_inertial  = R_b2i

        ephem = BinaryAsteroidEphemerides(times, r_sun, r_secondary, R_primary_to_secondary, R_primary_to_inertial)
        @test ephem isa BinaryAsteroidEphemerides{Vector{SMatrix{3,3,Float64,9}}}

        shape1  = load_shape_obj(path_obj; scale=1000, with_face_visibility=false, with_bvh=false)
        shape2  = load_shape_obj(path_obj; scale=500,  with_face_visibility=false, with_bvh=false)
        problem = BinaryAsteroidThermoPhysicalProblem(
            (shape1, shape2),
            (thermo_params, thermo_params);
            with_self_shadowing   = false,
            with_self_heating     = false,
            with_mutual_shadowing = false,
            with_mutual_heating   = false,
        )

        output = BinaryAsteroidOutputSpec(
            SingleAsteroidOutputSpec(output_times, subsurface_face_ids, true, true, false, true, true),
            SingleAsteroidOutputSpec(output_times, subsurface_face_ids, true, true, false, true, true),
        )
        solution = solve(problem, CrankNicolson();
            ephem                         = ephem,
            output                        = output,
            initial_temperature_primary   = 200.0,
            initial_temperature_secondary = 200.0,
            show_progress                 = true,
        )

        # Solution type
        @test solution isa BinaryAsteroidThermoPhysicalSolution
        @test solution.primary   isa SingleAsteroidThermoPhysicalSolution
        @test solution.secondary isa SingleAsteroidThermoPhysicalSolution

        # forces/torques are populated at output_times for both bodies
        @test length(solution.primary.forces)    == length(output_times)
        @test length(solution.secondary.torques) == length(output_times)

        # export: both bodies should have thermal_net_forces.csv with force/torque columns
        export_solution(DIR_OUTPUT, solution)
        for body in ("primary", "secondary")
            filepath = joinpath(DIR_OUTPUT, body, "thermal_net_forces.csv")
            df = CSV.read(filepath, DataFrame)
            @test "force_x"  in names(df)
            @test "torque_x" in names(df)
        end
    end
end
