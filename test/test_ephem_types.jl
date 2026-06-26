#=
test_ephem_types.jl

Unit tests for parametric *Ephemerides{R} types introduced in v0.2.0:
- Type hierarchy
- Construction via convenience constructors (plain Array input auto-converted)
- Automatic element-type conversion to SVector{3,Float64} / SMatrix{3,3,Float64,9}
- Explicit nothing constructor
- Field access
- Base.length
- Inner constructor dimension validation
- Explicit parametric construction {Nothing} and {<:AbstractVector}
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

    # Plain Julia arrays — constructors auto-convert to SVector/SMatrix internally.
    times                  = collect(range(0.0, 100.0; length=n))
    r_sun                  = [[au2m, 0.0, 0.0] for _ in 1:n]
    r_secondary            = [[1e6, 0.0, 0.0]  for _ in 1:n]
    R_body_to_inertial     = [Matrix{Float64}(I, 3, 3) for _ in 1:n]
    R_primary_to_secondary = [Matrix{Float64}(I, 3, 3) for _ in 1:n]
    R_primary_to_inertial  = [Matrix{Float64}(I, 3, 3) for _ in 1:n]

    @testset "Type hierarchy" begin
        ephem_s  = SingleAsteroidEphemerides(times, r_sun)
        ephem_sd = SingleAsteroidEphemerides(times, r_sun, R_body_to_inertial)
        ephem_b  = BinaryAsteroidEphemerides(times, r_sun, r_secondary, R_primary_to_secondary)
        ephem_bd = BinaryAsteroidEphemerides(times, r_sun, r_secondary, R_primary_to_secondary, R_primary_to_inertial)

        @test ephem_s  isa AbstractSingleAsteroidEphemerides
        @test ephem_sd isa AbstractSingleAsteroidEphemerides
        @test ephem_b  isa AbstractBinaryAsteroidEphemerides
        @test ephem_bd isa AbstractBinaryAsteroidEphemerides

        @test ephem_s  isa AbstractAsteroidEphemerides
        @test ephem_sd isa AbstractAsteroidEphemerides
        @test ephem_b  isa AbstractAsteroidEphemerides
        @test ephem_bd isa AbstractAsteroidEphemerides

        # Parametric type parameters
        @test ephem_s  isa SingleAsteroidEphemerides{Nothing}
        @test ephem_sd isa SingleAsteroidEphemerides{Vector{SMatrix{3,3,Float64,9}}}
        @test ephem_b  isa BinaryAsteroidEphemerides{Nothing}
        @test ephem_bd isa BinaryAsteroidEphemerides{Vector{SMatrix{3,3,Float64,9}}}
    end

    @testset "SingleAsteroidEphemerides{Nothing}" begin
        ephem = SingleAsteroidEphemerides(times, r_sun)

        @test ephem.times              === times
        @test ephem.r_sun              == r_sun
        @test eltype(ephem.r_sun)      === SVector{3, Float64}
        @test ephem.R_body_to_inertial === nothing
        @test length(ephem)            == n

        # DimensionMismatch
        @test_throws DimensionMismatch SingleAsteroidEphemerides{Nothing}(times[1:end-1], r_sun, nothing)
    end

    @testset "SingleAsteroidEphemerides{Nothing} — explicit nothing" begin
        ephem = SingleAsteroidEphemerides(times, r_sun, nothing)

        @test ephem isa SingleAsteroidEphemerides{Nothing}
        @test ephem.R_body_to_inertial === nothing
    end

    @testset "SingleAsteroidEphemerides — SMatrix fast path (no allocation)" begin
        # Passing Vector{SMatrix} directly hits the no-op fast path in _to_smat33_vec.
        R_smat = [SMatrix{3,3,Float64,9}(I) for _ in 1:n]
        ephem  = SingleAsteroidEphemerides(times, r_sun, R_smat)

        @test eltype(ephem.R_body_to_inertial)  === SMatrix{3, 3, Float64, 9}
        @test ephem.R_body_to_inertial          === R_smat  # same object: no copy
    end

    @testset "SingleAsteroidEphemerides with rotation matrices" begin
        ephem = SingleAsteroidEphemerides(times, r_sun, R_body_to_inertial)

        @test ephem.times                       === times
        @test ephem.r_sun                       == r_sun
        @test eltype(ephem.r_sun)               === SVector{3, Float64}
        @test ephem.R_body_to_inertial          == R_body_to_inertial
        @test eltype(ephem.R_body_to_inertial)  === SMatrix{3, 3, Float64, 9}
        @test length(ephem)                     == n

        # DimensionMismatch
        @test_throws DimensionMismatch SingleAsteroidEphemerides(times, r_sun[1:end-1], R_body_to_inertial)
        @test_throws DimensionMismatch SingleAsteroidEphemerides(times, r_sun, R_body_to_inertial[1:end-1])
    end

    @testset "BinaryAsteroidEphemerides{Nothing}" begin
        ephem = BinaryAsteroidEphemerides(times, r_sun, r_secondary, R_primary_to_secondary)

        @test ephem.times                          === times
        @test ephem.r_sun                          == r_sun
        @test eltype(ephem.r_sun)                  === SVector{3, Float64}
        @test ephem.r_secondary                    == r_secondary
        @test eltype(ephem.r_secondary)            === SVector{3, Float64}
        @test ephem.R_primary_to_secondary         == R_primary_to_secondary
        @test eltype(ephem.R_primary_to_secondary) === SMatrix{3, 3, Float64, 9}
        @test ephem.R_primary_to_inertial          === nothing
        @test length(ephem)                        == n

        # DimensionMismatch
        @test_throws DimensionMismatch BinaryAsteroidEphemerides{Nothing}(times[1:end-1], r_sun, r_secondary, R_primary_to_secondary, nothing)
        @test_throws DimensionMismatch BinaryAsteroidEphemerides(times, r_sun[1:end-1], r_secondary, R_primary_to_secondary)
        @test_throws DimensionMismatch BinaryAsteroidEphemerides(times, r_sun, r_secondary[1:end-1], R_primary_to_secondary)
        @test_throws DimensionMismatch BinaryAsteroidEphemerides(times, r_sun, r_secondary, R_primary_to_secondary[1:end-1])
    end

    @testset "BinaryAsteroidEphemerides{Nothing} — explicit nothing" begin
        ephem = BinaryAsteroidEphemerides(times, r_sun, r_secondary, R_primary_to_secondary, nothing)

        @test ephem isa BinaryAsteroidEphemerides{Nothing}
        @test ephem.R_primary_to_inertial === nothing
    end

    @testset "BinaryAsteroidEphemerides with rotation matrices" begin
        ephem = BinaryAsteroidEphemerides(times, r_sun, r_secondary, R_primary_to_secondary, R_primary_to_inertial)

        @test ephem.times                         === times
        @test ephem.r_sun                         == r_sun
        @test eltype(ephem.r_sun)                 === SVector{3, Float64}
        @test ephem.r_secondary                   == r_secondary
        @test eltype(ephem.r_secondary)           === SVector{3, Float64}
        @test ephem.R_primary_to_secondary        == R_primary_to_secondary
        @test eltype(ephem.R_primary_to_secondary) === SMatrix{3, 3, Float64, 9}
        @test ephem.R_primary_to_inertial         == R_primary_to_inertial
        @test eltype(ephem.R_primary_to_inertial) === SMatrix{3, 3, Float64, 9}
        @test length(ephem)                       == n

        # DimensionMismatch
        @test_throws DimensionMismatch BinaryAsteroidEphemerides(times, r_sun[1:end-1], r_secondary, R_primary_to_secondary, R_primary_to_inertial)
        @test_throws DimensionMismatch BinaryAsteroidEphemerides(times, r_sun, r_secondary[1:end-1], R_primary_to_secondary, R_primary_to_inertial)
        @test_throws DimensionMismatch BinaryAsteroidEphemerides(times, r_sun, r_secondary, R_primary_to_secondary[1:end-1], R_primary_to_inertial)
        @test_throws DimensionMismatch BinaryAsteroidEphemerides(times, r_sun, r_secondary, R_primary_to_secondary, R_primary_to_inertial[1:end-1])
    end

    @testset "AbstractRange as times" begin
        et_range = range(0.0, 100.0; length=n)

        ephem_s = SingleAsteroidEphemerides(et_range, r_sun)
        @test ephem_s.times isa Vector{Float64}
        @test ephem_s.times == collect(et_range)
        @test length(ephem_s) == n

        ephem_b = BinaryAsteroidEphemerides(et_range, r_sun, r_secondary, R_primary_to_secondary)
        @test ephem_b.times isa Vector{Float64}
        @test ephem_b.times == collect(et_range)
        @test length(ephem_b) == n
    end
end
