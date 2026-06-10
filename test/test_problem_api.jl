#=
test_problem_api.jl

Unit tests for the Problem-Solver API introduced in v0.2.0:
- subsolar_temperature
- BinaryAsteroidThermoPhysicalProblem convenience constructor (tuple args)
- init_temperature! with AbstractMatrix
- init_temperature! with per-body temperatures for binary systems
=#

@testset "Problem API" begin
    msg = """
    ┏━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━┓
    |                  Test: Problem API                     |
    ┗━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━┛
    """
    println(msg)

    shape1 = load_shape_obj(joinpath(@__DIR__, "shape", "single_face.obj"))
    shape2 = load_shape_obj(joinpath(@__DIR__, "shape", "icosahedron.obj"))
    thermo_params1 = ThermoParams(0.1, 1000.0, 700.0, 0.1, 0.0, 0.9, 0.1, 0.01, 10)
    thermo_params2 = ThermoParams(0.2, 1500.0, 700.0, 0.2, 0.0, 0.9, 0.1, 0.01, 10)

    @testset "subsolar_temperature" begin
        au2m = AsteroidThermoPhysicalModels.au2m
        r☉_1au  = SVector{3, Float64}(1 * au2m, 0, 0)
        r☉_2au  = SVector{3, Float64}(2 * au2m, 0, 0)

        Tss_1au = subsolar_temperature(r☉_1au, thermo_params1)
        Tss_2au = subsolar_temperature(r☉_2au, thermo_params1)

        @test Tss_1au > 0
        @test Tss_1au isa Float64
        @test subsolar_temperature(r☉_1au, 0.1, 0.9) ≈ Tss_1au  # two-argument form
        @test Tss_2au < Tss_1au                                 # farther from Sun → cooler
    end

    @testset "BinaryAsteroidThermoPhysicalProblem convenience constructor" begin
        problem = BinaryAsteroidThermoPhysicalProblem(
            (shape1, shape2),
            (thermo_params1, thermo_params2);
            with_self_shadowing   = false,
            with_self_heating     = false,
            with_mutual_shadowing = false,
            with_mutual_heating   = false,
        )

        @test problem isa BinaryAsteroidThermoPhysicalProblem
        @test problem.primary.shape   === shape1
        @test problem.secondary.shape === shape2

        # Single-body kwargs applied identically to both bodies
        @test problem.primary.with_self_shadowing == problem.secondary.with_self_shadowing == false
        @test problem.primary.with_self_heating   == problem.secondary.with_self_heating   == false
        @test problem.with_mutual_shadowing == false
        @test problem.with_mutual_heating   == false
        
        # Each body gets its own parameters
        @test problem.primary.thermo_params.thermal_conductivity[begin]   ≈ 0.1
        @test problem.secondary.thermo_params.thermal_conductivity[begin] ≈ 0.2
    end

    @testset "init_temperature! with AbstractMatrix" begin
        problem = SingleAsteroidThermoPhysicalProblem(shape1, thermo_params1;
            with_self_shadowing = false,
            with_self_heating   = false,
        )
        stpm = AsteroidThermoPhysicalModels._build_single_tpm(problem, CrankNicolson())

        n_depth  = thermo_params1.n_depth
        n_face   = length(shape1.faces)
        T_matrix = fill(350.0, n_depth, n_face)
        init_temperature!(stpm, T_matrix)

        @test stpm.temperature ≈ T_matrix
    end

    @testset "init_temperature! binary per-body" begin
        prob1 = SingleAsteroidThermoPhysicalProblem(shape1, thermo_params1;
            with_self_shadowing = false,
            with_self_heating   = false,
        )
        prob2 = SingleAsteroidThermoPhysicalProblem(shape2, thermo_params2;
            with_self_shadowing = false,
            with_self_heating   = false,
        )
        stpm1 = AsteroidThermoPhysicalModels._build_single_tpm(prob1, CrankNicolson())
        stpm2 = AsteroidThermoPhysicalModels._build_single_tpm(prob2, CrankNicolson())
        btpm  = BinaryAsteroidThermoPhysicalModel(stpm1, stpm2, false, false)

        init_temperature!(btpm, 200.0, 250.0)

        @test all(btpm.pri.temperature .== 200.0)
        @test all(btpm.sec.temperature .== 250.0)
    end
end
