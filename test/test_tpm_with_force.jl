#=
test_tpm_with_force.jl

Tests for the force/torque-enabled simulation paths via parametric types:
    SingleAsteroidEphemerides{<:AbstractVector}  →  SingleAsteroidThermoPhysicalSolution{Vector{SVector{3,Float64}}}
    BinaryAsteroidEphemerides{<:AbstractVector}  →  BinaryAsteroidThermoPhysicalSolution{Vector{SVector{3,Float64}}}

Covers: _alloc_solution_with_force, record_timestep! {<:AbstractVector},
        export_solution {<:AbstractVector}, _solve {<:AbstractVector} (single + binary).
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
    R_b2i = [SMatrix{3,3,Float64,9}(RotZ(2π * et / P)) for et in et_range]

    # Sun position in body-fixed frame (inverse rotation of inertial position)
    r_sun = [SVector{3,Float64}(inv(RotZ(2π * et / P)) * SVector(au2m, 0.0, 0.0)) for et in et_range]

    path_obj      = joinpath("shape", "icosahedron.obj")
    thermo_params = ThermoParams(0.1, 1000.0, 600.0, 0.1, 0.0, 0.9, 0.1, 0.01, 5)
    times_to_save = times[end-n_step_cycle:end]
    face_ID       = [1]

    @testset "Single asteroid — {<:AbstractVector}" begin
        DIR_OUTPUT = mktempdir()

        ephem = SingleAsteroidEphemerides(times, r_sun, R_b2i)
        @test ephem isa SingleAsteroidEphemerides{Vector{SMatrix{3,3,Float64,9}}}

        shape   = load_shape_obj(path_obj; scale=1000, with_face_visibility=false, with_bvh=false)
        problem = AsteroidThermoPhysicalModels.SingleAsteroidThermoPhysicalProblem(shape, thermo_params;
            with_self_shadowing = false,
            with_self_heating   = false,
        )

        output = SingleAsteroidOutputSpec(times_to_save, face_ID)
        solution = solve(problem, CrankNicolson();
            ephem               = ephem,
            output              = output,
            initial_temperature = 200.0,
            show_progress       = true,
        )

        # Solution type
        @test solution isa SingleAsteroidThermoPhysicalSolution{Vector{SVector{3,Float64}}}

        # forces/torques are populated
        @test solution.forces  isa Vector{SVector{3,Float64}}
        @test solution.torques isa Vector{SVector{3,Float64}}
        @test length(solution.forces)  == length(times)
        @test length(solution.torques) == length(times)

        # export: physical_quantities.csv must contain force/torque columns
        AsteroidThermoPhysicalModels.export_solution(DIR_OUTPUT, solution)
        @test isfile(joinpath(DIR_OUTPUT, "physical_quantities.csv"))
        df = CSV.read(joinpath(DIR_OUTPUT, "physical_quantities.csv"), DataFrame)
        @test "force_x"  in names(df)
        @test "force_y"  in names(df)
        @test "force_z"  in names(df)
        @test "torque_x" in names(df)
        @test "torque_y" in names(df)
        @test "torque_z" in names(df)
    end

    @testset "Binary asteroid — {<:AbstractVector}" begin
        DIR_OUTPUT = mktempdir()

        r_secondary            = [SVector{3,Float64}(1e4, 0, 0) for _ in et_range]
        R_primary_to_secondary = [SMatrix{3,3,Float64,9}(I)     for _ in et_range]
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
            SingleAsteroidOutputSpec(times_to_save, face_ID),
            SingleAsteroidOutputSpec(times_to_save, face_ID),
        )
        solution = solve(problem, CrankNicolson();
            ephem                         = ephem,
            output                        = output,
            initial_temperature_primary   = 200.0,
            initial_temperature_secondary = 200.0,
            show_progress                 = true,
        )

        # Solution type
        @test solution isa BinaryAsteroidThermoPhysicalSolution{Vector{SVector{3,Float64}}}
        @test solution.primary   isa SingleAsteroidThermoPhysicalSolution{Vector{SVector{3,Float64}}}
        @test solution.secondary isa SingleAsteroidThermoPhysicalSolution{Vector{SVector{3,Float64}}}

        # forces/torques are populated for both bodies
        @test length(solution.primary.forces)    == length(times)
        @test length(solution.secondary.torques) == length(times)

        # export: both bodies should have force/torque columns
        AsteroidThermoPhysicalModels.export_solution(DIR_OUTPUT, solution)
        for body in ("primary", "secondary")
            df = CSV.read(joinpath(DIR_OUTPUT, body, "physical_quantities.csv"), DataFrame)
            @test "force_x"  in names(df)
            @test "torque_x" in names(df)
        end
    end
end
