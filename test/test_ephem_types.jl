#=
test_ephem_types.jl

Unit tests for *Ephemerides types introduced in v0.2.0:
- Type hierarchy
- Construction and field access
- Base.length
- Inner constructor dimension validation
- Upgrade constructors
=#

@testset "Ephemerides types" begin
    msg = """
    ┏━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━┓
    |               Test: Ephemerides types                  |
    ┗━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━┛
    """
    println(msg)

    n = 10
    au2m = AsteroidThermoPhysicalModels.au2m

    times                  = collect(range(0.0, 100.0; length=n))
    r_sun                  = [SVector{3, Float64}(au2m, 0, 0) for _ in 1:n]
    r_secondary            = [SVector{3, Float64}(1e6, 0, 0)  for _ in 1:n]
    R_body_to_inertial     = [SMatrix{3,3,Float64,9}(I) for _ in 1:n]
    R_primary_to_secondary = [SMatrix{3,3,Float64,9}(I) for _ in 1:n]
    R_primary_to_inertial  = [SMatrix{3,3,Float64,9}(I) for _ in 1:n]

    @testset "Type hierarchy" begin
        ephem_s  = SingleAsteroidEphemerides(times, r_sun)
        ephem_sd = SingleAsteroidEphemeridesWithDynamics(times, r_sun, R_body_to_inertial)
        ephem_b  = BinaryAsteroidEphemerides(times, r_sun, r_secondary, R_primary_to_secondary)
        ephem_bd = BinaryAsteroidEphemeridesWithDynamics(times, r_sun, r_secondary, R_primary_to_secondary, R_primary_to_inertial)

        @test ephem_s  isa AbstractSingleAsteroidEphemerides
        @test ephem_sd isa AbstractSingleAsteroidEphemerides
        @test ephem_b  isa AbstractBinaryAsteroidEphemerides
        @test ephem_bd isa AbstractBinaryAsteroidEphemerides

        @test ephem_s  isa AbstractAsteroidEphemerides
        @test ephem_sd isa AbstractAsteroidEphemerides
        @test ephem_b  isa AbstractAsteroidEphemerides
        @test ephem_bd isa AbstractAsteroidEphemerides
    end

    @testset "SingleAsteroidEphemerides" begin
        ephem = SingleAsteroidEphemerides(times, r_sun)

        @test ephem.times === times
        @test ephem.r_sun === r_sun
        @test length(ephem) == n

        # DimensionMismatch
        @test_throws DimensionMismatch SingleAsteroidEphemerides(times[1:end-1], r_sun)
    end

    @testset "SingleAsteroidEphemeridesWithDynamics" begin
        ephem = SingleAsteroidEphemeridesWithDynamics(times, r_sun, R_body_to_inertial)

        @test ephem.times              === times
        @test ephem.r_sun              === r_sun
        @test ephem.R_body_to_inertial === R_body_to_inertial
        @test length(ephem) == n

        @test_throws DimensionMismatch SingleAsteroidEphemeridesWithDynamics(times, r_sun[1:end-1], R_body_to_inertial)
        @test_throws DimensionMismatch SingleAsteroidEphemeridesWithDynamics(times, r_sun, R_body_to_inertial[1:end-1])

        # Upgrade constructor
        ephem_base = SingleAsteroidEphemerides(times, r_sun)
        ephem_up   = SingleAsteroidEphemeridesWithDynamics(ephem_base, R_body_to_inertial)

        @test ephem_up isa SingleAsteroidEphemeridesWithDynamics
        @test ephem_up.times === ephem_base.times
        @test ephem_up.r_sun === ephem_base.r_sun
        @test ephem_up.R_body_to_inertial === R_body_to_inertial
    end

    @testset "BinaryAsteroidEphemerides" begin
        ephem = BinaryAsteroidEphemerides(times, r_sun, r_secondary, R_primary_to_secondary)

        @test ephem.times                  === times
        @test ephem.r_sun                  === r_sun
        @test ephem.r_secondary            === r_secondary
        @test ephem.R_primary_to_secondary === R_primary_to_secondary
        @test length(ephem) == n

        @test_throws DimensionMismatch BinaryAsteroidEphemerides(times, r_sun[1:end-1], r_secondary, R_primary_to_secondary)
        @test_throws DimensionMismatch BinaryAsteroidEphemerides(times, r_sun, r_secondary[1:end-1], R_primary_to_secondary)
        @test_throws DimensionMismatch BinaryAsteroidEphemerides(times, r_sun, r_secondary, R_primary_to_secondary[1:end-1])
    end

    @testset "BinaryAsteroidEphemeridesWithDynamics" begin
        ephem = BinaryAsteroidEphemeridesWithDynamics(times, r_sun, r_secondary, R_primary_to_secondary, R_primary_to_inertial)

        @test ephem.times                  === times
        @test ephem.r_sun                  === r_sun
        @test ephem.r_secondary            === r_secondary
        @test ephem.R_primary_to_secondary === R_primary_to_secondary
        @test ephem.R_primary_to_inertial  === R_primary_to_inertial
        @test length(ephem) == n

        @test_throws DimensionMismatch BinaryAsteroidEphemeridesWithDynamics(times, r_sun[1:end-1], r_secondary, R_primary_to_secondary, R_primary_to_inertial)
        @test_throws DimensionMismatch BinaryAsteroidEphemeridesWithDynamics(times, r_sun, r_secondary[1:end-1], R_primary_to_secondary, R_primary_to_inertial)
        @test_throws DimensionMismatch BinaryAsteroidEphemeridesWithDynamics(times, r_sun, r_secondary, R_primary_to_secondary[1:end-1], R_primary_to_inertial)
        @test_throws DimensionMismatch BinaryAsteroidEphemeridesWithDynamics(times, r_sun, r_secondary, R_primary_to_secondary, R_primary_to_inertial[1:end-1])

        # Upgrade constructor
        ephem_base = BinaryAsteroidEphemerides(times, r_sun, r_secondary, R_primary_to_secondary)
        ephem_up   = BinaryAsteroidEphemeridesWithDynamics(ephem_base, R_primary_to_inertial)

        @test ephem_up isa BinaryAsteroidEphemeridesWithDynamics
        @test ephem_up.times                  === ephem_base.times
        @test ephem_up.r_sun                  === ephem_base.r_sun
        @test ephem_up.r_secondary            === ephem_base.r_secondary
        @test ephem_up.R_primary_to_secondary === ephem_base.R_primary_to_secondary
        @test ephem_up.R_primary_to_inertial  === R_primary_to_inertial
    end
end
