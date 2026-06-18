#=
test_output_spec.jl

Unit tests for SingleAsteroidOutputSpec and BinaryAsteroidOutputSpec types
introduced in v0.2.0:
- Construction and field access
- _validate_output_spec: valid case (no error)
- _validate_output_spec: invalid case (ArgumentError when times_to_save ⊄ ephem.times)
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

    times_to_save = times[end-2:end]
    face_ID       = [1, 3, 5]

    @testset "SingleAsteroidOutputSpec construction" begin
        output = SingleAsteroidOutputSpec(times_to_save, face_ID)

        @test output isa SingleAsteroidOutputSpec
        @test output.times_to_save === times_to_save
        @test output.face_ID       === face_ID
    end

    @testset "BinaryAsteroidOutputSpec construction" begin
        spec1 = SingleAsteroidOutputSpec(times_to_save, face_ID)
        spec2 = SingleAsteroidOutputSpec(times_to_save, [2, 4])
        output   = BinaryAsteroidOutputSpec(spec1, spec2)

        @test output isa BinaryAsteroidOutputSpec
        @test output.primary   === spec1
        @test output.secondary === spec2
    end

    @testset "_validate_output_spec — valid (no error)" begin
        output = SingleAsteroidOutputSpec(times_to_save, face_ID)
        @test_nowarn AsteroidThermoPhysicalModels._validate_output_spec(output, ephem)
    end

    @testset "_validate_output_spec — invalid single (ArgumentError)" begin
        invalid_time  = times[end] + 1.0  # not in ephem.times
        output_bad    = SingleAsteroidOutputSpec([times[1], invalid_time], face_ID)
        @test_throws ArgumentError AsteroidThermoPhysicalModels._validate_output_spec(output_bad, ephem)
    end

    @testset "_validate_output_spec — valid binary (no error)" begin
        r_secondary            = [SVector{3, Float64}(1e6, 0, 0) for _ in 1:n]
        R_primary_to_secondary = [SMatrix{3,3,Float64,9}(I)      for _ in 1:n]
        ephem_b = BinaryAsteroidEphemerides(times, r_sun, r_secondary, R_primary_to_secondary)

        output = BinaryAsteroidOutputSpec(
            SingleAsteroidOutputSpec(times_to_save, face_ID),
            SingleAsteroidOutputSpec(times_to_save, [2]),
        )
        @test_nowarn AsteroidThermoPhysicalModels._validate_output_spec(output, ephem_b)
    end

    @testset "_validate_output_spec — invalid binary secondary (ArgumentError)" begin
        r_secondary            = [SVector{3, Float64}(1e6, 0, 0) for _ in 1:n]
        R_primary_to_secondary = [SMatrix{3,3,Float64,9}(I)      for _ in 1:n]
        ephem_b = BinaryAsteroidEphemerides(times, r_sun, r_secondary, R_primary_to_secondary)

        invalid_time = times[end] + 1.0
        output = BinaryAsteroidOutputSpec(
            SingleAsteroidOutputSpec(times_to_save, face_ID),
            SingleAsteroidOutputSpec([invalid_time], [1]),
        )
        @test_throws ArgumentError AsteroidThermoPhysicalModels._validate_output_spec(output, ephem_b)
    end
end