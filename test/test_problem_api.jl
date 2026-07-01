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

    thermo_params1 = ThermoParams(
        conductivity    = 0.1,
        density         = 1000.0,
        heat_capacity   = 700.0,
        reflectance_vis = 0.1,
        reflectance_ir  = 0.0,
        emissivity      = 0.9,
    )
    thermo_params2 = ThermoParams(
        conductivity    = 0.2,
        density         = 1500.0,
        heat_capacity   = 700.0,
        reflectance_vis = 0.2,
        reflectance_ir  = 0.0,
        emissivity      = 0.9,
    )
    grid_params1 = GridParams(; z_max=0.1, n_depth=10)
    grid_params2 = GridParams(; z_max=0.1, n_depth=10)

    @testset "subsolar_temperature" begin
        au2m = AsteroidThermoPhysicalModels.au2m
        r☉_1au  = [1 * au2m, 0.0, 0.0]
        r☉_2au  = [2 * au2m, 0.0, 0.0]

        Tss_1au = subsolar_temperature(r☉_1au, 0.1, 0.9)
        Tss_2au = subsolar_temperature(r☉_2au, 0.1, 0.9)

        @test Tss_1au > 0
        @test Tss_1au isa Float64
        @test Tss_2au < Tss_1au  # farther from Sun → cooler
    end

    @testset "BinaryAsteroidThermoPhysicalProblem convenience constructor" begin
        problem = BinaryAsteroidThermoPhysicalProblem(
            (shape1, shape2),
            (thermo_params1, thermo_params2),
            (grid_params1, grid_params2);
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
        @test problem.primary.thermo_params.conductivity[begin]   ≈ 0.1
        @test problem.secondary.thermo_params.conductivity[begin] ≈ 0.2
    end

    @testset "init_temperature! with AbstractMatrix" begin
        problem = SingleAsteroidThermoPhysicalProblem(shape1, thermo_params1, grid_params1;
            with_self_shadowing = false,
            with_self_heating   = false,
        )
        state = AsteroidThermoPhysicalModels._build_single_state(problem, CrankNicolson())

        n_depth  = grid_params1.n_depth
        n_face   = length(shape1.faces)
        T_matrix = fill(350.0, n_depth, n_face)
        AsteroidThermoPhysicalModels.init_temperature!(state, T_matrix)

        @test state.temperature ≈ T_matrix
    end

    @testset "_expand_thermo_params — wrong vector length" begin
        # icosahedron has 20 faces; a ThermoParams with 3 entries is neither 1 nor n_face
        tp_wrong = ThermoParams([0.1, 0.2, 0.3], [1000.0, 1000.0, 1000.0], [700.0, 700.0, 700.0],
                                [0.1, 0.1, 0.1], [0.0, 0.0, 0.0], [0.9, 0.9, 0.9])
        @test_throws ArgumentError SingleAsteroidThermoPhysicalProblem(shape2, tp_wrong, grid_params2;
            with_self_shadowing = false,
            with_self_heating   = false,
        )
    end

    @testset "init_temperature! binary per-body" begin
        prob = BinaryAsteroidThermoPhysicalProblem(
            (shape1, shape2),
            (thermo_params1, thermo_params2),
            (grid_params1, grid_params2);
            with_self_shadowing   = false,
            with_self_heating     = false,
            with_mutual_shadowing = false,
            with_mutual_heating   = false,
        )
        state = AsteroidThermoPhysicalModels._build_binary_state(prob, CrankNicolson())

        AsteroidThermoPhysicalModels.init_temperature!(state, 200.0, 250.0)

        @test all(state.primary.temperature   .== 200.0)
        @test all(state.secondary.temperature .== 250.0)

        AsteroidThermoPhysicalModels.init_temperature!(state, 300.0)

        @test all(state.primary.temperature   .== 300.0)
        @test all(state.secondary.temperature .== 300.0)
    end
end
