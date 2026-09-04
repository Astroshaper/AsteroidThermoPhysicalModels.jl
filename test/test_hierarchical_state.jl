#=
test_hierarchical_state.jl

Unit tests for HierarchicalSingleAsteroidThermoPhysicalState:
- _build_single_state dispatch for HierarchicalShapeModel
- init_temperature! for HierarchicalSingleAsteroidThermoPhysicalState (Real and AbstractMatrix)
- automatic preparation of the geometric data required for self-shadowing and self-heating
- update_flux_sun! on the global level, compared against the plain ShapeModel result
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

    @testset "self-shadowing geometry prepared automatically" begin
        # Self-shadowing is evaluated on the global shape, so the visibility graph and the
        # maximum elevations must be built there rather than on the hierarchical wrapper.
        hier_shape3 = load_shape_obj(joinpath(@__DIR__, "shape", "icosahedron.obj"); as_hierarchical=true)
        add_roughness_models!(hier_shape3, roughness_model)

        @test isnothing(hier_shape3.global_shape.face_visibility_graph)
        @test isnothing(hier_shape3.global_shape.face_max_elevations)

        SingleAsteroidThermoPhysicalProblem(hier_shape3, thermo_params, grid_params;
            with_self_shadowing = true,
            with_self_heating   = false,
        )

        @test !isnothing(hier_shape3.global_shape.face_visibility_graph)
        @test !isnothing(hier_shape3.global_shape.face_max_elevations)
    end

    @testset "self-heating geometry prepared automatically" begin
        hier_shape4 = load_shape_obj(joinpath(@__DIR__, "shape", "icosahedron.obj"); as_hierarchical=true)
        add_roughness_models!(hier_shape4, roughness_model)

        @test isnothing(hier_shape4.global_shape.face_visibility_graph)

        SingleAsteroidThermoPhysicalProblem(hier_shape4, thermo_params, grid_params;
            with_self_shadowing = false,
            with_self_heating   = true,
        )

        @test !isnothing(hier_shape4.global_shape.face_visibility_graph)
        @test isnothing(hier_shape4.global_shape.face_max_elevations)  # not needed by self-heating
    end

    @testset "update_flux_sun! (global level)" begin
        # The global level of a hierarchical state must reproduce the plain ShapeModel result
        # exactly: the sub-face machinery must not perturb the global solar flux.
        path_obj = joinpath(@__DIR__, "shape", "icosahedron.obj")
        au2m = AsteroidThermoPhysicalModels.au2m
        r☉s = [
            SVector(1.0, 0.0, 0.0) * au2m,        # 1.0 au
            SVector(0.3, -0.5, 0.8) * 1.2au2m,    # 1.2 au, oblique
        ]

        for with_self_shadowing in (false, true)
            shape_plain = load_shape_obj(path_obj)
            shape_hier  = load_shape_obj(path_obj; as_hierarchical=true)
            add_roughness_models!(shape_hier, roughness_model)

            problem_plain = SingleAsteroidThermoPhysicalProblem(shape_plain, thermo_params, grid_params;
                with_self_shadowing, with_self_heating=false)
            problem_hier  = SingleAsteroidThermoPhysicalProblem(shape_hier, thermo_params, grid_params;
                with_self_shadowing, with_self_heating=false)

            state_plain = AsteroidThermoPhysicalModels._build_single_state(problem_plain, CrankNicolson())
            state_hier  = AsteroidThermoPhysicalModels._build_single_state(problem_hier,  CrankNicolson())

            for r☉ in r☉s
                AsteroidThermoPhysicalModels.update_flux_sun!(state_plain, r☉)
                AsteroidThermoPhysicalModels.update_flux_sun!(state_hier,  r☉)

                # Global level matches the plain ShapeModel result exactly
                @test state_hier.illuminated_faces == state_plain.illuminated_faces
                @test state_hier.flux_sun          == state_plain.flux_sun

                # Physical sanity: the icosahedron is convex, so illumination reduces to n̂ ⋅ r̂☉ > 0.
                # Faces whose normal is perpendicular to the Sun direction (cosθ ≈ 0 up to
                # round-off) are skipped, since their illumination flag is not well defined.
                r̂☉ = normalize(r☉)
                F☉ = AsteroidThermoPhysicalModels.SOLAR_CONST / (norm(r☉) * AsteroidThermoPhysicalModels.m2au)^2
                for (i, n̂) in enumerate(shape_hier.global_shape.face_normals)
                    cosθ = n̂ ⋅ r̂☉
                    abs(cosθ) < 1e-8 && continue
                    if cosθ > 0
                        @test state_hier.illuminated_faces[i]
                        @test state_hier.flux_sun[i] ≈ F☉ * cosθ
                    else
                        @test !state_hier.illuminated_faces[i]
                        @test state_hier.flux_sun[i] == 0.0
                    end
                end
                @test any(state_hier.illuminated_faces) && !all(state_hier.illuminated_faces)

                # Sub-face states are not touched by the global-level update
                for rs in state_hier.roughness_states
                    @test !any(rs.illuminated_faces)
                    @test all(rs.flux_sun .== 0.0)
                end
            end
        end
    end

    @testset "update_flux_sun! requires face_visibility_graph for self-shadowing" begin
        # The problem constructor builds the graph, but the global shape is mutable and may
        # lose it afterwards. Self-shadowing must then fail loudly rather than silently skip
        # the shadow test.
        shape_hier = load_shape_obj(joinpath(@__DIR__, "shape", "icosahedron.obj"); as_hierarchical=true)
        add_roughness_models!(shape_hier, roughness_model)
        problem_hier = SingleAsteroidThermoPhysicalProblem(shape_hier, thermo_params, grid_params;
            with_self_shadowing=true, with_self_heating=false)
        state_hier = AsteroidThermoPhysicalModels._build_single_state(problem_hier, CrankNicolson())

        shape_hier.global_shape.face_visibility_graph = nothing
        r☉ = SVector(1.0, 0.0, 0.0) * AsteroidThermoPhysicalModels.au2m
        @test_throws ErrorException AsteroidThermoPhysicalModels.update_flux_sun!(state_hier, r☉)
    end
end
