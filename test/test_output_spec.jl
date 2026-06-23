#=
test_output_spec.jl

Unit tests for SingleAsteroidOutputSpec and BinaryAsteroidOutputSpec types
introduced in v0.2.0:
- Construction and field access (including Bool flag defaults)
- Inner constructor validation (ArgumentError when save_subsurface_temperature=true with empty subsurface_face_ids)
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
    r_sun  = [SVector{3, Float64}(au2m, 0, 0) for _ in 1:n]
    ephem  = SingleAsteroidEphemerides(times, r_sun)

    output_times        = times[end-2:end]
    subsurface_face_ids = [1, 3, 5]

    @testset "SingleAsteroidOutputSpec construction" begin
        output = SingleAsteroidOutputSpec(output_times, subsurface_face_ids)

        @test output isa SingleAsteroidOutputSpec
        @test output.output_times        === output_times
        @test output.subsurface_face_ids === subsurface_face_ids

        # Default Bool flags
        @test output.save_surface_temperature    == true
        @test output.save_subsurface_temperature == true
        @test output.save_face_forces            == false
        @test output.save_forces                 == true
        @test output.save_torques                == true
    end

    @testset "SingleAsteroidOutputSpec inner constructor validation" begin
        # save_subsurface_temperature=true with empty subsurface_face_ids → ArgumentError
        @test_throws ArgumentError SingleAsteroidOutputSpec(output_times, Int[], true, true, false, true, true)
        # save_subsurface_temperature=false with empty subsurface_face_ids → OK
        @test_nowarn SingleAsteroidOutputSpec(output_times, Int[], true, false, false, true, true)
    end

    @testset "BinaryAsteroidOutputSpec construction" begin
        spec1 = SingleAsteroidOutputSpec(output_times, subsurface_face_ids)
        spec2 = SingleAsteroidOutputSpec(output_times, [2, 4])
        output   = BinaryAsteroidOutputSpec(spec1, spec2)

        @test output isa BinaryAsteroidOutputSpec
        @test output.primary   === spec1
        @test output.secondary === spec2
    end

    @testset "_validate_output_spec — valid (no error)" begin
        output = SingleAsteroidOutputSpec(output_times, subsurface_face_ids)
        @test_nowarn AsteroidThermoPhysicalModels._validate_output_spec(output, ephem)
    end

    @testset "_validate_output_spec — invalid single (ArgumentError)" begin
        invalid_time  = times[end] + 1.0  # not in ephem.times
        output_bad    = SingleAsteroidOutputSpec([times[1], invalid_time], subsurface_face_ids)
        @test_throws ArgumentError AsteroidThermoPhysicalModels._validate_output_spec(output_bad, ephem)
    end

    @testset "_validate_output_spec — valid binary (no error)" begin
        r_secondary            = [SVector{3, Float64}(1e6, 0, 0) for _ in 1:n]
        R_primary_to_secondary = [SMatrix{3,3,Float64,9}(I)      for _ in 1:n]
        ephem_b = BinaryAsteroidEphemerides(times, r_sun, r_secondary, R_primary_to_secondary)

        output = BinaryAsteroidOutputSpec(
            SingleAsteroidOutputSpec(output_times, subsurface_face_ids),
            SingleAsteroidOutputSpec(output_times, [2]),
        )
        @test_nowarn AsteroidThermoPhysicalModels._validate_output_spec(output, ephem_b)
    end

    @testset "_validate_output_spec — invalid binary secondary (ArgumentError)" begin
        r_secondary            = [SVector{3, Float64}(1e6, 0, 0) for _ in 1:n]
        R_primary_to_secondary = [SMatrix{3,3,Float64,9}(I)      for _ in 1:n]
        ephem_b = BinaryAsteroidEphemerides(times, r_sun, r_secondary, R_primary_to_secondary)

        invalid_time = times[end] + 1.0
        output = BinaryAsteroidOutputSpec(
            SingleAsteroidOutputSpec(output_times, subsurface_face_ids),
            SingleAsteroidOutputSpec([invalid_time], [1]),
        )
        @test_throws ArgumentError AsteroidThermoPhysicalModels._validate_output_spec(output, ephem_b)
    end
end