#=
test_hierarchical_state.jl

Unit tests for HierarchicalSingleAsteroidThermoPhysicalState:
- _build_single_state dispatch for HierarchicalShapeModel
- init_temperature! for HierarchicalSingleAsteroidThermoPhysicalState (Real and AbstractMatrix)
=#

@testset "HierarchicalSingleAsteroidThermoPhysicalState" begin
    msg = """
    ┏━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━┓
    |        Test: HierarchicalSingleAsteroidTPMState        |
    ┗━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━┛
    """
    println(msg)

    # Build a HierarchicalShapeModel from icosahedron with a crater roughness on all faces
    hier_shape = load_shape_obj(joinpath(@__DIR__, "shape", "icosahedron.obj"); as_hierarchical=true)
    roughness_model = create_shape_crater(0.4, 0.1; Nx=4, Ny=4)
    add_roughness_models!(hier_shape, roughness_model)

    thermo_params = ThermoParams(
        conductivity    = 0.1,
        density         = 1000.0,
        heat_capacity   = 700.0,
        reflectance_vis = 0.1,
        reflectance_ir  = 0.0,
        emissivity      = 0.9,
    )
    grid_params = GridParams(; z_max=0.1, n_depth=10)

    problem = SingleAsteroidThermoPhysicalProblem(hier_shape, thermo_params, grid_params;
        with_self_shadowing = false,
        with_self_heating   = false,
    )

    n_global_faces = length(hier_shape.global_shape.faces)

    @testset "problem construction" begin
        @test problem isa SingleAsteroidThermoPhysicalProblem
        @test problem.shape === hier_shape
        @test length(problem.thermo_params.conductivity) == n_global_faces
    end

    state = AsteroidThermoPhysicalModels._build_single_state(problem, CrankNicolson())

    @testset "_build_single_state returns HierarchicalSingleAsteroidThermoPhysicalState" begin
        @test state isa AsteroidThermoPhysicalModels.HierarchicalSingleAsteroidThermoPhysicalState
        @test size(state.temperature) == (grid_params.n_depth, n_global_faces)
        @test length(state.illuminated_faces) == n_global_faces
        @test length(state.flux_sun)          == n_global_faces
        @test length(state.flux_scat)         == n_global_faces
        @test length(state.flux_rad)          == n_global_faces
        @test length(state.face_forces)       == n_global_faces
    end

    @testset "face_roughness_indices and roughness_states" begin
        @test length(state.face_roughness_indices) == n_global_faces
        # All faces have roughness → every index is positive and sequential
        @test all(state.face_roughness_indices .> 0)
        @test state.face_roughness_indices == 1:n_global_faces
        @test length(state.roughness_states) == n_global_faces
        # Each sub-state is an independent SingleAsteroidThermoPhysicalState
        @test all(rs isa AsteroidThermoPhysicalModels.SingleAsteroidThermoPhysicalState
                  for rs in state.roughness_states)
    end

    @testset "init_temperature! (Real)" begin
        AsteroidThermoPhysicalModels.init_temperature!(state, 200.0)

        @test all(state.temperature .== 200.0)
        @test all(state.roughness_states[1].temperature .== 200.0)
        @test all(state.roughness_states[end].temperature .== 200.0)
    end

    @testset "init_temperature! (AbstractMatrix)" begin
        T_mat = fill(300.0, grid_params.n_depth, n_global_faces)
        AsteroidThermoPhysicalModels.init_temperature!(state, T_mat)

        @test state.temperature ≈ T_mat
        # Sub-faces initialized to surface temperature of parent global face (300.0)
        @test all(state.roughness_states[1].temperature .== 300.0)
        @test all(state.roughness_states[end].temperature .== 300.0)
    end

    @testset "surface_temperature" begin
        AsteroidThermoPhysicalModels.init_temperature!(state, 250.0)
        T_surf = AsteroidThermoPhysicalModels.surface_temperature(state)

        @test T_surf isa Vector{Float64}
        @test length(T_surf) == n_global_faces
        @test all(T_surf .== 250.0)
    end

    @testset "no-roughness faces" begin
        # Build a shape with roughness on only one face
        hier_shape2 = load_shape_obj(joinpath(@__DIR__, "shape", "icosahedron.obj"); as_hierarchical=true)
        add_roughness_models!(hier_shape2, roughness_model, 1)

        problem2 = SingleAsteroidThermoPhysicalProblem(hier_shape2, thermo_params, grid_params;
            with_self_shadowing = false,
            with_self_heating   = false,
        )
        state2 = AsteroidThermoPhysicalModels._build_single_state(problem2, CrankNicolson())

        n_faces2 = length(hier_shape2.global_shape.faces)

        @test length(state2.face_roughness_indices) == n_faces2
        @test state2.face_roughness_indices[1] == 1       # face 1 has roughness
        @test all(state2.face_roughness_indices[2:end] .== 0)  # others don't
        @test length(state2.roughness_states) == 1
    end
end
