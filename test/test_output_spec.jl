#=
test_output_spec.jl

Unit tests for SingleAsteroidOutputSpec and BinaryAsteroidOutputSpec types:
- Keyword construction and defaults; face-selected outputs are switched on by listing faces
- Input conversion (ranges and integer vectors of other element types)
- Binary convenience constructor
- _validate_output_spec: valid case (no error)
- _validate_output_spec: invalid case (ArgumentError when output_times ⊄ ephem.times)
=#

@testset "OutputSpec types" begin
    msg = """
    ┏━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━┓
    |               Test: OutputSpec types                   |
    ┗━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━┛
    """
    println(msg)

    au2m = AsteroidThermoPhysicalModels.au2m

    n = 10
    times  = collect(range(0.0, 100.0; length=n))
    r_sun  = [[au2m, 0.0, 0.0] for _ in 1:n]
    ephem  = SingleAsteroidEphemerides(times, r_sun)

    output_times        = times[end-2:end]
    subsurface_face_ids = [1, 3, 5]

    @testset "SingleAsteroidOutputSpec defaults" begin
        output = SingleAsteroidOutputSpec(output_times)

        @test output isa SingleAsteroidOutputSpec
        @test output.output_times == output_times

        # Only the surface temperature is saved by default
        @test output.save_surface_temperature == true
        @test output.subsurface_face_ids      == Int[]
        @test output.roughness_face_ids       == Int[]
        @test output.save_face_forces         == false
        @test output.save_forces              == false
        @test output.save_torques             == false
    end

    @testset "SingleAsteroidOutputSpec keywords" begin
        output = SingleAsteroidOutputSpec(output_times;
            save_surface_temperature = false,
            subsurface_face_ids      = subsurface_face_ids,
            roughness_face_ids       = [1, 7],
            save_face_forces         = true,
            save_forces              = true,
            save_torques             = true,
        )

        @test output.save_surface_temperature == false
        @test output.subsurface_face_ids      == subsurface_face_ids
        @test output.roughness_face_ids       == [1, 7]
        @test output.save_face_forces         == true
        @test output.save_forces              == true
        @test output.save_torques             == true
    end

    @testset "SingleAsteroidOutputSpec input conversion" begin
        # Ranges and other integer element types are converted to the field types
        output = SingleAsteroidOutputSpec(range(0.0, 100.0; length=n);
            subsurface_face_ids = 1:3, roughness_face_ids = Int32[2, 4])

        @test output.output_times        isa Vector{Float64}
        @test output.subsurface_face_ids isa Vector{Int}
        @test output.roughness_face_ids  isa Vector{Int}
        @test output.subsurface_face_ids == [1, 2, 3]
        @test output.roughness_face_ids  == [2, 4]
    end

    @testset "BinaryAsteroidOutputSpec construction" begin
        spec1 = SingleAsteroidOutputSpec(output_times; subsurface_face_ids)
        spec2 = SingleAsteroidOutputSpec(output_times; subsurface_face_ids=[2, 4])
        output   = BinaryAsteroidOutputSpec(spec1, spec2)

        @test output isa BinaryAsteroidOutputSpec
        @test output.primary   === spec1
        @test output.secondary === spec2
    end

    @testset "_validate_output_spec — valid (no error)" begin
        output = SingleAsteroidOutputSpec(output_times; subsurface_face_ids)
        @test_nowarn AsteroidThermoPhysicalModels._validate_output_spec(output, ephem)
    end

    @testset "_validate_output_spec — invalid single (ArgumentError: time not in ephem)" begin
        invalid_time  = times[end] + 1.0  # not in ephem.times
        output_bad    = SingleAsteroidOutputSpec([times[1], invalid_time]; subsurface_face_ids)
        @test_throws ArgumentError AsteroidThermoPhysicalModels._validate_output_spec(output_bad, ephem)
    end

    @testset "_validate_output_spec — invalid single (ArgumentError: forces with no rotation matrices)" begin
        # save_forces=true requires R_body_to_inertial in ephemerides
        output_with_forces = SingleAsteroidOutputSpec(output_times; subsurface_face_ids, save_forces=true)
        @test_throws ArgumentError AsteroidThermoPhysicalModels._validate_output_spec(output_with_forces, ephem)
    end

    @testset "_validate_output_spec — valid binary (no error)" begin
        r_secondary            = [[1e6, 0.0, 0.0]          for _ in 1:n]
        R_primary_to_secondary = [Matrix{Float64}(I, 3, 3) for _ in 1:n]
        ephem_b = BinaryAsteroidEphemerides(times, r_sun, r_secondary, R_primary_to_secondary)

        output = BinaryAsteroidOutputSpec(
            SingleAsteroidOutputSpec(output_times; subsurface_face_ids),
            SingleAsteroidOutputSpec(output_times; subsurface_face_ids=[2]),
        )
        @test_nowarn AsteroidThermoPhysicalModels._validate_output_spec(output, ephem_b)
    end

    @testset "_validate_output_spec — invalid binary secondary (ArgumentError)" begin
        r_secondary            = [[1e6, 0.0, 0.0]          for _ in 1:n]
        R_primary_to_secondary = [Matrix{Float64}(I, 3, 3) for _ in 1:n]
        ephem_b = BinaryAsteroidEphemerides(times, r_sun, r_secondary, R_primary_to_secondary)

        invalid_time = times[end] + 1.0
        output = BinaryAsteroidOutputSpec(
            SingleAsteroidOutputSpec(output_times; subsurface_face_ids),
            SingleAsteroidOutputSpec([invalid_time]; subsurface_face_ids=[1]),
        )
        @test_throws ArgumentError AsteroidThermoPhysicalModels._validate_output_spec(output, ephem_b)
    end

    @testset "BinaryAsteroidOutputSpec convenience constructor" begin
        output_times_s        = times[end-1:end]
        subsurface_face_ids_s = [2, 4]

        output = BinaryAsteroidOutputSpec(output_times, output_times_s;
            subsurface_face_ids_primary   = subsurface_face_ids,
            subsurface_face_ids_secondary = subsurface_face_ids_s,
            save_surface_temperature = false,
            save_face_forces         = true,
            save_forces              = false,
            save_torques             = false,
        )

        @test output isa BinaryAsteroidOutputSpec

        # Each body gets the correct output_times and subsurface_face_ids
        @test output.primary.output_times          == output_times
        @test output.primary.subsurface_face_ids   == subsurface_face_ids
        @test output.secondary.output_times        == output_times_s
        @test output.secondary.subsurface_face_ids == subsurface_face_ids_s

        # Shared Bool flags are propagated identically to both bodies
        @test output.primary.save_surface_temperature   == false
        @test output.secondary.save_surface_temperature == false
        @test output.primary.save_face_forces   == true
        @test output.secondary.save_face_forces == true
        @test output.primary.save_forces        == false
        @test output.secondary.save_forces      == false
        @test output.primary.save_torques       == false
        @test output.secondary.save_torques     == false
    end
end
