#=
test_hierarchical_state.jl

Unit tests for HierarchicalSingleAsteroidThermoPhysicalState:
- _build_single_state dispatch for HierarchicalShapeModel
- init_temperature! for HierarchicalSingleAsteroidThermoPhysicalState (Real and AbstractMatrix)
- automatic preparation of the geometric data required for self-shadowing and self-heating
- update_flux_sun! on the global level, compared against the plain ShapeModel result
- update_flux_sun! on the sub-face level: gate on the global illumination, consistency with a
  standalone roughness model, and a near-flat roughness model reproducing the parent face
- update_flux_scat_single! / update_flux_rad_single! on the global level, likewise
- update_flux_scat_single! / update_flux_rad_single! on the sub-face level: self-heating inside
  the roughness model, the external irradiation from the other global faces, and its absence
  when the global self-heating is off
- update_temperature! on the global level, for all three solvers and for zero conductivity
- update_temperature! on the sub-face level: the sub-faces advance, a near-flat roughness model
  reproduces the temperature of its parent face, and zero conductivity gives radiative
  equilibrium on every sub-face
- update_thermal_force!: faces without roughness keep the plain result, a near-flat roughness
  model reproduces the force on its parent face, the roughness scale cancels, a deep crater
  changes the force and the net force/torque are re-summed, and the sky re-absorption term
  follows the global self-heating flag
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

                # Sub-faces of a dark global face are all dark (the sub-face level is covered
                # in detail by the dedicated testsets below)
                for (i, k) in enumerate(state_hier.face_roughness_indices)
                    k == 0 && continue
                    state_hier.illuminated_faces[i] && continue
                    rs = state_hier.roughness_states[k]
                    @test !any(rs.illuminated_faces)
                    @test all(rs.flux_sun .== 0.0)
                end
            end
        end
    end

    @testset "update_flux_sun! (sub-face level)" begin
        # For a lit global face, the sub-faces must be exactly what a standalone state of the
        # roughness model gives for the Sun vector rotated into the local frame. For a dark
        # global face, every sub-face must be dark regardless of the local geometry: the
        # roughness model has no terrain beyond its rim, so a Sun just below the local horizon
        # would otherwise light sub-faces tilted towards it.
        path_obj = joinpath(@__DIR__, "shape", "icosahedron.obj")
        crater   = create_shape_crater(0.4, 0.1; Nx=6, Ny=6)

        shape_hier = load_shape_obj(path_obj; as_hierarchical=true)
        add_roughness_models!(shape_hier, crater)
        problem_hier = SingleAsteroidThermoPhysicalProblem(shape_hier, thermo_params, grid_params;
            with_self_shadowing=true, with_self_heating=false)
        state_hier = AsteroidThermoPhysicalModels._build_single_state(problem_hier, CrankNicolson())

        # A standalone state of the same roughness model, illuminated with the same model
        problem_alone = SingleAsteroidThermoPhysicalProblem(crater, thermo_params, grid_params;
            with_self_shadowing=true, with_self_heating=false)
        state_alone = AsteroidThermoPhysicalModels._build_single_state(problem_alone, CrankNicolson())

        # Sub-face self-shadowing is always on, so the roughness model carries its geometry
        @test !isnothing(crater.face_visibility_graph)
        @test !isnothing(crater.face_max_elevations)

        r☉s = [
            SVector(1.0, 0.0, 0.0) * AsteroidThermoPhysicalModels.au2m,
            SVector(0.3, -0.5, 0.8) * 1.2AsteroidThermoPhysicalModels.au2m,
        ]
        for r☉ in r☉s
            AsteroidThermoPhysicalModels.update_flux_sun!(state_hier, r☉)

            n_lit_parents = n_dark_parents = 0
            for (i, k) in enumerate(state_hier.face_roughness_indices)
                rs = state_hier.roughness_states[k]
                if state_hier.illuminated_faces[i]
                    n_lit_parents += 1
                    r☉_local = transform_physical_vector_global_to_local(shape_hier, i, r☉)
                    @test norm(r☉_local) ≈ norm(r☉)   # pure rotation: solar flux is preserved
                    AsteroidThermoPhysicalModels.update_flux_sun!(state_alone, r☉_local)
                    @test rs.illuminated_faces == state_alone.illuminated_faces
                    @test rs.flux_sun          == state_alone.flux_sun
                    @test any(rs.illuminated_faces)
                else
                    n_dark_parents += 1
                    @test !any(rs.illuminated_faces)
                    @test all(rs.flux_sun .== 0.0)
                end
            end
            # Both branches of the gate were exercised
            @test n_lit_parents > 0 && n_dark_parents > 0
        end
    end

    @testset "update_flux_sun! (near-flat roughness reproduces the parent face)" begin
        # With a roughness model that is flat to 1e-6, every sub-face has the normal of its
        # parent face, so its solar flux must equal the parent's. This pins the rotation into
        # the local frame: a wrong rotation would change cos θ and show up here.
        shape_hier = load_shape_obj(joinpath(@__DIR__, "shape", "icosahedron.obj"); as_hierarchical=true)
        add_roughness_models!(shape_hier, create_shape_crater(0.4, 1e-6; Nx=4, Ny=4))
        problem_hier = SingleAsteroidThermoPhysicalProblem(shape_hier, thermo_params, grid_params;
            with_self_shadowing=false, with_self_heating=false)
        state_hier = AsteroidThermoPhysicalModels._build_single_state(problem_hier, CrankNicolson())

        r☉ = SVector(0.3, -0.2, 0.9) * AsteroidThermoPhysicalModels.au2m
        AsteroidThermoPhysicalModels.update_flux_sun!(state_hier, r☉)

        # The residual tilt of the sub-face normals (~1e-6) is an absolute error in cos θ, so
        # compare against the flux at normal incidence rather than relatively: at grazing
        # incidence the parent flux itself is tiny and a relative tolerance is ill-conditioned.
        # A wrong rotation would produce O(1) differences and is still caught.
        r̂☉ = normalize(r☉)
        F☉ = AsteroidThermoPhysicalModels.SOLAR_CONST / (norm(r☉) * AsteroidThermoPhysicalModels.m2au)^2
        n_checked = 0
        for (i, k) in enumerate(state_hier.face_roughness_indices)
            rs = state_hier.roughness_states[k]
            if state_hier.illuminated_faces[i]
                shape_hier.global_shape.face_normals[i] ⋅ r̂☉ < 1e-3 && continue  # grazing: skip
                n_checked += 1
                @test all(rs.illuminated_faces)
                @test all(isapprox.(rs.flux_sun, state_hier.flux_sun[i]; atol=1e-5 * F☉))
            else
                @test all(rs.flux_sun .== 0.0)
            end
        end
        @test n_checked > 0
    end

    @testset "update_flux_sun! (partial roughness)" begin
        # Faces without a roughness model are skipped; the one with it is updated
        shape_hier = load_shape_obj(joinpath(@__DIR__, "shape", "icosahedron.obj"); as_hierarchical=true)
        add_roughness_models!(shape_hier, create_shape_crater(0.4, 0.1; Nx=4, Ny=4), 1)
        problem_hier = SingleAsteroidThermoPhysicalProblem(shape_hier, thermo_params, grid_params;
            with_self_shadowing=false, with_self_heating=false)
        state_hier = AsteroidThermoPhysicalModels._build_single_state(problem_hier, CrankNicolson())

        # Point the Sun along the normal of face 1 so that it is certainly lit
        n̂₁ = shape_hier.global_shape.face_normals[1]
        AsteroidThermoPhysicalModels.update_flux_sun!(state_hier, n̂₁ * AsteroidThermoPhysicalModels.au2m)

        @test state_hier.illuminated_faces[1]
        @test length(state_hier.roughness_states) == 1
        @test all(state_hier.roughness_states[1].illuminated_faces)
        @test all(>(0), state_hier.roughness_states[1].flux_sun)
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

    @testset "update_flux_scat_single! / update_flux_rad_single! (global level)" begin
        # Self-heating vanishes identically on a convex shape, so the icosahedron used above
        # cannot tell a working implementation from a broken one. A crater is concave: its
        # faces see each other, and both self-heating terms are non-zero.
        make_crater(; as_hierarchical) =
            create_shape_crater(0.4, 0.1; Nx=8, Ny=8, as_hierarchical)

        shape_plain = make_crater(as_hierarchical=false)
        shape_hier  = make_crater(as_hierarchical=true)
        add_roughness_models!(shape_hier, create_shape_crater(0.4, 0.1; Nx=4, Ny=4))

        problem_plain = SingleAsteroidThermoPhysicalProblem(shape_plain, thermo_params, grid_params;
            with_self_shadowing=false, with_self_heating=true)
        problem_hier  = SingleAsteroidThermoPhysicalProblem(shape_hier, thermo_params, grid_params;
            with_self_shadowing=false, with_self_heating=true)

        state_plain = AsteroidThermoPhysicalModels._build_single_state(problem_plain, CrankNicolson())
        state_hier  = AsteroidThermoPhysicalModels._build_single_state(problem_hier,  CrankNicolson())

        # Vary the temperature from face to face, so that a wrong face index in the radiation
        # term would change the result rather than cancel out.
        n_faces = length(shape_plain.faces)
        T₀ = repeat(reshape(range(200.0, 300.0; length=n_faces) |> collect, 1, n_faces),
                    grid_params.n_depth)
        AsteroidThermoPhysicalModels.init_temperature!(state_plain, T₀)
        AsteroidThermoPhysicalModels.init_temperature!(state_hier,  T₀)

        r☉ = SVector(0.2, -0.1, 1.0) * AsteroidThermoPhysicalModels.au2m
        for state in (state_plain, state_hier)
            AsteroidThermoPhysicalModels.update_flux_sun!(state, r☉)
            AsteroidThermoPhysicalModels.update_flux_scat_single!(state)
            AsteroidThermoPhysicalModels.update_flux_rad_single!(state)
        end

        # Global level matches the plain ShapeModel result exactly
        @test state_hier.flux_scat == state_plain.flux_scat
        @test state_hier.flux_rad  == state_plain.flux_rad

        # The comparison is not vacuous: the concave shape really does exchange energy
        @test any(>(0), state_hier.flux_scat)
        @test any(>(0), state_hier.flux_rad)

        # Sub-faces are updated as well (covered in detail by the dedicated testsets below)
        for rs in state_hier.roughness_states
            @test any(>(0), rs.flux_rad)
        end
    end

    @testset "update_flux_scat_single! / update_flux_rad_single! (sub-face level)" begin
        # Self-heating inside the roughness model must be exactly what a standalone state of
        # the roughness model gives, and the flux the parent receives from the other global
        # faces must be added on top, weighted by each sub-face's sky view factor. A concave
        # global crater makes the parent fluxes non-zero; a roughness crater makes the
        # intra-model exchange non-zero.
        crater = create_shape_crater(0.4, 0.1; Nx=6, Ny=6)
        shape_hier = create_shape_crater(0.4, 0.1; Nx=8, Ny=8, as_hierarchical=true)
        add_roughness_models!(shape_hier, crater)
        problem_hier = SingleAsteroidThermoPhysicalProblem(shape_hier, thermo_params, grid_params;
            with_self_shadowing=false, with_self_heating=true)
        state_hier = AsteroidThermoPhysicalModels._build_single_state(problem_hier, CrankNicolson())

        # Self-heating inside the roughness model is always on
        @test all(rs.problem.with_self_heating for rs in state_hier.roughness_states)
        @test !isnothing(crater.face_visibility_graph)

        # Standalone state of the same roughness model, with the same flags as the sub-states
        problem_alone = SingleAsteroidThermoPhysicalProblem(crater, thermo_params, grid_params;
            with_self_shadowing=true, with_self_heating=true)
        state_alone = AsteroidThermoPhysicalModels._build_single_state(problem_alone, CrankNicolson())

        # Vary the temperature from face to face so the radiation term is not degenerate
        n_faces = length(shape_hier.global_shape.faces)
        T₀ = repeat(reshape(range(200.0, 300.0; length=n_faces) |> collect, 1, n_faces), grid_params.n_depth)
        AsteroidThermoPhysicalModels.init_temperature!(state_hier, T₀)

        r☉ = SVector(0.2, -0.1, 1.0) * AsteroidThermoPhysicalModels.au2m
        AsteroidThermoPhysicalModels.update_flux_sun!(state_hier, r☉)
        AsteroidThermoPhysicalModels.update_flux_scat_single!(state_hier)
        AsteroidThermoPhysicalModels.update_flux_rad_single!(state_hier)

        # The parent fluxes are non-zero somewhere, so the external term is exercised
        @test any(>(0), state_hier.flux_scat)
        @test any(>(0), state_hier.flux_rad)

        # Sky view factor of each sub-face of the (shared) roughness model
        f_sky = [max(0.0, 1 - sum(get_view_factors(crater.face_visibility_graph, j)))
                 for j in eachindex(crater.faces)]
        @test all(0 .<= f_sky .<= 1)
        @test any(<(1), f_sky)   # the crater walls do hide part of the sky

        n_checked = 0
        for (i, k) in enumerate(state_hier.face_roughness_indices)
            rs = state_hier.roughness_states[k]
            state_hier.illuminated_faces[i] || continue
            n_checked += 1

            # Standalone: same illumination and temperature, self-heating inside the model only
            r☉_local = transform_physical_vector_global_to_local(shape_hier, i, r☉)
            AsteroidThermoPhysicalModels.init_temperature!(state_alone, T₀[begin, i])
            AsteroidThermoPhysicalModels.update_flux_sun!(state_alone, r☉_local)
            AsteroidThermoPhysicalModels.update_flux_scat_single!(state_alone)
            AsteroidThermoPhysicalModels.update_flux_rad_single!(state_alone)

            @test rs.flux_scat ≈ state_alone.flux_scat .+ f_sky .* state_hier.flux_scat[i]
            @test rs.flux_rad  ≈ state_alone.flux_rad  .+ f_sky .* state_hier.flux_rad[i]
        end
        @test n_checked > 0
    end

    @testset "global self-heating off removes the external term only" begin
        # With `with_self_heating = false` on the global problem, the parent fluxes are zero
        # and the sub-faces receive nothing from outside — but the self-heating inside each
        # roughness model is still computed. This is what makes an on/off comparison of the
        # global self-heating isolate the external contribution.
        crater = create_shape_crater(0.4, 0.1; Nx=6, Ny=6)
        shape_hier = create_shape_crater(0.4, 0.1; Nx=8, Ny=8, as_hierarchical=true)
        add_roughness_models!(shape_hier, crater)
        problem_hier = SingleAsteroidThermoPhysicalProblem(shape_hier, thermo_params, grid_params;
            with_self_shadowing=false, with_self_heating=false)
        state_hier = AsteroidThermoPhysicalModels._build_single_state(problem_hier, CrankNicolson())
        AsteroidThermoPhysicalModels.init_temperature!(state_hier, 250.0)

        problem_alone = SingleAsteroidThermoPhysicalProblem(crater, thermo_params, grid_params;
            with_self_shadowing=true, with_self_heating=true)
        state_alone = AsteroidThermoPhysicalModels._build_single_state(problem_alone, CrankNicolson())
        AsteroidThermoPhysicalModels.init_temperature!(state_alone, 250.0)

        r☉ = SVector(0.2, -0.1, 1.0) * AsteroidThermoPhysicalModels.au2m
        AsteroidThermoPhysicalModels.update_flux_sun!(state_hier, r☉)
        AsteroidThermoPhysicalModels.update_flux_scat_single!(state_hier)
        AsteroidThermoPhysicalModels.update_flux_rad_single!(state_hier)

        @test all(state_hier.flux_scat .== 0.0)
        @test all(state_hier.flux_rad  .== 0.0)

        for (i, k) in enumerate(state_hier.face_roughness_indices)
            rs = state_hier.roughness_states[k]
            state_hier.illuminated_faces[i] || continue
            AsteroidThermoPhysicalModels.update_flux_sun!(state_alone, transform_physical_vector_global_to_local(shape_hier, i, r☉))
            AsteroidThermoPhysicalModels.update_flux_scat_single!(state_alone)
            AsteroidThermoPhysicalModels.update_flux_rad_single!(state_alone)
            @test rs.flux_scat == state_alone.flux_scat
            @test rs.flux_rad  == state_alone.flux_rad
            @test any(>(0), rs.flux_rad)   # the crater does heat itself
        end
    end

    @testset "near-flat roughness receives exactly the parent flux from outside" begin
        # A flat roughness model has no walls (sky view factor 1) and no exchange between its
        # own faces, so each sub-face must end up with exactly the parent's flux.
        shape_hier = create_shape_crater(0.4, 0.1; Nx=8, Ny=8, as_hierarchical=true)
        add_roughness_models!(shape_hier, create_shape_crater(0.4, 1e-6; Nx=4, Ny=4))
        problem_hier = SingleAsteroidThermoPhysicalProblem(shape_hier, thermo_params, grid_params;
            with_self_shadowing=false, with_self_heating=true)
        state_hier = AsteroidThermoPhysicalModels._build_single_state(problem_hier, CrankNicolson())
        AsteroidThermoPhysicalModels.init_temperature!(state_hier, 250.0)

        r☉ = SVector(0.2, -0.1, 1.0) * AsteroidThermoPhysicalModels.au2m
        AsteroidThermoPhysicalModels.update_flux_sun!(state_hier, r☉)
        AsteroidThermoPhysicalModels.update_flux_scat_single!(state_hier)
        AsteroidThermoPhysicalModels.update_flux_rad_single!(state_hier)

        # The view factors of a model flat to 1e-6 are ~1e-11 rather than zero, so the exchange
        # between its own faces leaves a residual of ~1e-10 W/m². Compare with an absolute
        # tolerance scaled by the largest parent flux: a relative one would fail on parents
        # whose flux is exactly zero.
        atol_scat = 1e-6 * maximum(state_hier.flux_scat)
        atol_rad  = 1e-6 * maximum(state_hier.flux_rad)
        @test atol_scat > 0 && atol_rad > 0
        for (i, k) in enumerate(state_hier.face_roughness_indices)
            rs = state_hier.roughness_states[k]
            @test all(isapprox.(rs.flux_scat, state_hier.flux_scat[i]; atol=atol_scat))
            @test all(isapprox.(rs.flux_rad,  state_hier.flux_rad[i];  atol=atol_rad))
        end
    end

    @testset "self-heating disabled leaves the global fluxes at zero" begin
        shape_hier = create_shape_crater(0.4, 0.1; Nx=8, Ny=8, as_hierarchical=true)
        add_roughness_models!(shape_hier, create_shape_crater(0.4, 0.1; Nx=4, Ny=4))

        problem_hier = SingleAsteroidThermoPhysicalProblem(shape_hier, thermo_params, grid_params;
            with_self_shadowing=false, with_self_heating=false)
        state_hier = AsteroidThermoPhysicalModels._build_single_state(problem_hier, CrankNicolson())
        AsteroidThermoPhysicalModels.init_temperature!(state_hier, 300.0)

        AsteroidThermoPhysicalModels.update_flux_sun!(state_hier, SVector(0.0, 0.0, 1.0) * AsteroidThermoPhysicalModels.au2m)
        AsteroidThermoPhysicalModels.update_flux_scat_single!(state_hier)
        AsteroidThermoPhysicalModels.update_flux_rad_single!(state_hier)

        @test all(state_hier.flux_scat .== 0.0)
        @test all(state_hier.flux_rad  .== 0.0)
        @test any(>(0), state_hier.flux_sun)  # the shape is lit, so the zeros are the flag's doing
    end

    @testset "update_temperature! (global level)" begin
        # The global faces of a hierarchical state are solved independently of their roughness
        # models, so after any number of steps they must match the plain ShapeModel exactly.
        # A concave crater with self-heating exercises every flux term in the surface balance.
        make_crater(; as_hierarchical) =
            create_shape_crater(0.4, 0.1; Nx=8, Ny=8, as_hierarchical)

        n_faces = length(make_crater(as_hierarchical=false).faces)
        T₀ = repeat(reshape(range(200.0, 300.0; length=n_faces) |> collect, 1, n_faces),
                    grid_params.n_depth)
        r☉ = SVector(0.2, -0.1, 1.0) * AsteroidThermoPhysicalModels.au2m
        Δt = 100.0  # λ = αΔt/Δz² ≈ 0.12 keeps the explicit Euler step stable
        n_steps = 10

        for algorithm in (CrankNicolson(), ImplicitEuler(), ExplicitEuler())
            shape_plain = make_crater(as_hierarchical=false)
            shape_hier  = make_crater(as_hierarchical=true)
            add_roughness_models!(shape_hier, create_shape_crater(0.4, 0.1; Nx=4, Ny=4))

            problem_plain = SingleAsteroidThermoPhysicalProblem(shape_plain, thermo_params, grid_params;
                with_self_shadowing=false, with_self_heating=true)
            problem_hier  = SingleAsteroidThermoPhysicalProblem(shape_hier, thermo_params, grid_params;
                with_self_shadowing=false, with_self_heating=true)

            state_plain = AsteroidThermoPhysicalModels._build_single_state(problem_plain, algorithm)
            state_hier  = AsteroidThermoPhysicalModels._build_single_state(problem_hier,  algorithm)
            AsteroidThermoPhysicalModels.init_temperature!(state_plain, T₀)
            AsteroidThermoPhysicalModels.init_temperature!(state_hier,  T₀)

            for _ in 1:n_steps, state in (state_plain, state_hier)
                AsteroidThermoPhysicalModels.update_flux_sun!(state, r☉)
                AsteroidThermoPhysicalModels.update_flux_scat_single!(state)
                AsteroidThermoPhysicalModels.update_flux_rad_single!(state)
                AsteroidThermoPhysicalModels.update_temperature!(state, Δt)
            end

            # Global level matches the plain ShapeModel result exactly, and has actually moved
            @test state_hier.temperature == state_plain.temperature
            @test state_hier.temperature != T₀

            # Sub-faces advance too (every face radiates, so none stays at its initial value)
            for (i, k) in enumerate(state_hier.face_roughness_indices)
                k == 0 && continue
                @test any(state_hier.roughness_states[k].temperature .!= T₀[begin, i])
            end
        end
    end

    @testset "update_temperature! (sub-face level, near-flat reproduces the parent)" begin
        # A roughness model flat to 1e-6 sees the same flux as its parent face and has the same
        # material and grid, so after any number of steps every one of its columns must match
        # the parent's column. This exercises the whole chain — local illumination, sub-face
        # boundary condition, per-state solver cache — against the global-level result, which
        # #226 pinned to the plain ShapeModel. The convex icosahedron keeps the external
        # irradiation at zero, so nothing else enters.
        path_obj = joinpath(@__DIR__, "shape", "icosahedron.obj")
        r☉ = SVector(0.3, -0.2, 0.9) * AsteroidThermoPhysicalModels.au2m
        r̂☉ = normalize(r☉)
        Δt = 100.0
        n_steps = 10

        for algorithm in (CrankNicolson(), ImplicitEuler(), ExplicitEuler())
            shape_hier = load_shape_obj(path_obj; as_hierarchical=true)
            add_roughness_models!(shape_hier, create_shape_crater(0.4, 1e-6; Nx=4, Ny=4))
            problem_hier = SingleAsteroidThermoPhysicalProblem(shape_hier, thermo_params, grid_params;
                with_self_shadowing=true, with_self_heating=false)
            state_hier = AsteroidThermoPhysicalModels._build_single_state(problem_hier, algorithm)
            AsteroidThermoPhysicalModels.init_temperature!(state_hier, 250.0)

            for _ in 1:n_steps
                AsteroidThermoPhysicalModels.update_flux_sun!(state_hier, r☉)
                AsteroidThermoPhysicalModels.update_flux_scat_single!(state_hier)
                AsteroidThermoPhysicalModels.update_flux_rad_single!(state_hier)
                AsteroidThermoPhysicalModels.update_temperature!(state_hier, Δt)
            end
            @test any(state_hier.temperature .!= 250.0)

            n_checked = 0
            for (i, k) in enumerate(state_hier.face_roughness_indices)
                rs = state_hier.roughness_states[k]
                # Skip grazing incidence, where the 1e-6 tilt of the sub-face normals is not
                # small against cos θ of the parent (same reasoning as for the fluxes)
                state_hier.illuminated_faces[i] && shape_hier.global_shape.face_normals[i] ⋅ r̂☉ < 0.05 && continue
                n_checked += 1
                for j in axes(rs.temperature, 2)
                    @test all(isapprox.(rs.temperature[:, j], state_hier.temperature[:, i]; rtol=1e-5))
                end
            end
            @test n_checked > 0
        end
    end

    @testset "update_temperature! (global level, zero conductivity)" begin
        # Zero conductivity replaces the solver by instantaneous radiative equilibrium at the
        # surface, which is the one code path that used to reach for `shape.faces` directly.
        thermo_params_k0 = ThermoParams(
            conductivity    = 0.0,
            density         = 1000.0,
            heat_capacity   = 700.0,
            reflectance_vis = 0.1,
            reflectance_ir  = 0.0,
            emissivity      = 0.9,
        )

        shape_plain = create_shape_crater(0.4, 0.1; Nx=8, Ny=8)
        shape_hier  = create_shape_crater(0.4, 0.1; Nx=8, Ny=8, as_hierarchical=true)
        add_roughness_models!(shape_hier, create_shape_crater(0.4, 0.1; Nx=4, Ny=4))

        problem_plain = SingleAsteroidThermoPhysicalProblem(shape_plain, thermo_params_k0, grid_params;
            with_self_shadowing=false, with_self_heating=true)
        problem_hier  = SingleAsteroidThermoPhysicalProblem(shape_hier, thermo_params_k0, grid_params;
            with_self_shadowing=false, with_self_heating=true)

        state_plain = AsteroidThermoPhysicalModels._build_single_state(problem_plain, CrankNicolson())
        state_hier  = AsteroidThermoPhysicalModels._build_single_state(problem_hier,  CrankNicolson())
        AsteroidThermoPhysicalModels.init_temperature!(state_plain, 250.0)
        AsteroidThermoPhysicalModels.init_temperature!(state_hier,  250.0)

        r☉ = SVector(0.2, -0.1, 1.0) * AsteroidThermoPhysicalModels.au2m
        for state in (state_plain, state_hier)
            AsteroidThermoPhysicalModels.update_flux_sun!(state, r☉)
            AsteroidThermoPhysicalModels.update_flux_scat_single!(state)
            AsteroidThermoPhysicalModels.update_flux_rad_single!(state)
            AsteroidThermoPhysicalModels.update_temperature!(state, 100.0)
        end

        @test state_hier.temperature == state_plain.temperature

        # Surface is in radiative equilibrium with the absorbed flux: εσT⁴ = F_abs
        εσ = 0.9 * AsteroidThermoPhysicalModels.σ_SB
        for i in axes(state_hier.temperature, 2)
            F_abs = AsteroidThermoPhysicalModels.absorbed_energy_flux(
                0.1, 0.0, state_hier.flux_sun[i], state_hier.flux_scat[i], state_hier.flux_rad[i])
            @test state_hier.temperature[begin, i] ≈ (F_abs / εσ)^(1/4)
        end
        @test any(>(0), state_hier.temperature[begin, :])

        # The same holds on every sub-face, with the sub-face's own fluxes
        for rs in state_hier.roughness_states
            for j in axes(rs.temperature, 2)
                F_abs = AsteroidThermoPhysicalModels.absorbed_energy_flux(
                    0.1, 0.0, rs.flux_sun[j], rs.flux_scat[j], rs.flux_rad[j])
                @test rs.temperature[begin, j] ≈ (F_abs / εσ)^(1/4)
            end
            @test any(>(0), rs.temperature[begin, :])
        end
    end

    # ---- update_thermal_force! ----------------------------------------------------------

    # One flux update followed by the thermal force, on a state whose temperature is set.
    function update_fluxes_and_force!(state, r☉)
        AsteroidThermoPhysicalModels.update_flux_sun!(state, r☉)
        AsteroidThermoPhysicalModels.update_flux_scat_single!(state)
        AsteroidThermoPhysicalModels.update_flux_rad_single!(state)
        AsteroidThermoPhysicalModels.update_thermal_force!(state)
    end

    @testset "update_thermal_force! (global level, partial roughness)" begin
        # A global face without a roughness model keeps the plain ShapeModel force exactly;
        # the one face with a crater gets a different force (its recoil is that of the crater).
        path_obj = joinpath(@__DIR__, "shape", "icosahedron.obj")
        shape_plain = load_shape_obj(path_obj)
        shape_hier  = load_shape_obj(path_obj; as_hierarchical=true)
        add_roughness_models!(shape_hier, roughness_model, 1)

        problem_plain = SingleAsteroidThermoPhysicalProblem(shape_plain, thermo_params, grid_params;
            with_self_shadowing=false, with_self_heating=false)
        problem_hier  = SingleAsteroidThermoPhysicalProblem(shape_hier, thermo_params, grid_params;
            with_self_shadowing=false, with_self_heating=false)
        state_plain = AsteroidThermoPhysicalModels._build_single_state(problem_plain, CrankNicolson())
        state_hier  = AsteroidThermoPhysicalModels._build_single_state(problem_hier,  CrankNicolson())
        AsteroidThermoPhysicalModels.init_temperature!(state_plain, 250.0)
        AsteroidThermoPhysicalModels.init_temperature!(state_hier,  250.0)

        r☉ = SVector(0.3, -0.2, 0.9) * AsteroidThermoPhysicalModels.au2m
        update_fluxes_and_force!(state_plain, r☉)
        update_fluxes_and_force!(state_hier,  r☉)

        @test state_hier.face_forces[2:end] == state_plain.face_forces[2:end]
        @test state_hier.face_forces[1]     != state_plain.face_forces[1]
        @test state_hier.force ≈ sum(state_hier.face_forces)
    end

    @testset "update_thermal_force! (near-flat roughness reproduces the parent face)" begin
        # A roughness model flat to 1e-6 has no self-heating and the normal of its parent, so
        # the sum of its sub-face forces, counted for the parent's area and rotated back, must
        # be the plain force on the parent. This pins both the rotation and the A_i / A_proj
        # normalisation.
        path_obj = joinpath(@__DIR__, "shape", "icosahedron.obj")
        shape_plain = load_shape_obj(path_obj)
        shape_hier  = load_shape_obj(path_obj; as_hierarchical=true)
        add_roughness_models!(shape_hier, create_shape_crater(0.4, 1e-6; Nx=4, Ny=4))

        problem_plain = SingleAsteroidThermoPhysicalProblem(shape_plain, thermo_params, grid_params;
            with_self_shadowing=false, with_self_heating=false)
        problem_hier  = SingleAsteroidThermoPhysicalProblem(shape_hier, thermo_params, grid_params;
            with_self_shadowing=false, with_self_heating=false)
        state_plain = AsteroidThermoPhysicalModels._build_single_state(problem_plain, CrankNicolson())
        state_hier  = AsteroidThermoPhysicalModels._build_single_state(problem_hier,  CrankNicolson())
        AsteroidThermoPhysicalModels.init_temperature!(state_plain, 250.0)
        AsteroidThermoPhysicalModels.init_temperature!(state_hier,  250.0)

        r☉ = SVector(0.3, -0.2, 0.9) * AsteroidThermoPhysicalModels.au2m
        r̂☉ = normalize(r☉)
        update_fluxes_and_force!(state_plain, r☉)
        update_fluxes_and_force!(state_hier,  r☉)

        n_checked = 0
        for i in eachindex(shape_plain.faces)
            # Skip grazing incidence, where the 1e-6 tilt of the sub-face normals is not small
            # against cos θ of the parent (same reasoning as for the fluxes)
            state_hier.illuminated_faces[i] && shape_plain.face_normals[i] ⋅ r̂☉ < 0.05 && continue
            n_checked += 1
            @test isapprox(state_hier.face_forces[i], state_plain.face_forces[i]; rtol=1e-4)
        end
        @test n_checked > 0
        @test isapprox(state_hier.force, state_plain.force; rtol=1e-4)
    end

    @testset "update_thermal_force! (roughness scale cancels)" begin
        # The force on one patch grows as scale², the number of patches covering the face falls
        # as scale⁻²: the face force must not depend on the scale of the roughness model.
        path_obj = joinpath(@__DIR__, "shape", "icosahedron.obj")
        r☉ = SVector(0.3, -0.2, 0.9) * AsteroidThermoPhysicalModels.au2m

        states = map((1.0, 0.1)) do scale
            shape_hier = load_shape_obj(path_obj; as_hierarchical=true)
            add_roughness_models!(shape_hier, roughness_model; scale)
            problem_hier = SingleAsteroidThermoPhysicalProblem(shape_hier, thermo_params, grid_params;
                with_self_shadowing=false, with_self_heating=false)
            state_hier = AsteroidThermoPhysicalModels._build_single_state(problem_hier, CrankNicolson())
            AsteroidThermoPhysicalModels.init_temperature!(state_hier, 250.0)
            update_fluxes_and_force!(state_hier, r☉)
            state_hier
        end

        @test all(isapprox.(states[1].face_forces, states[2].face_forces; rtol=1e-12))
        @test isapprox(states[1].force,  states[2].force;  rtol=1e-12)
        @test isapprox(states[1].torque, states[2].torque; rtol=1e-12)
        @test any(f -> norm(f) > 0, states[1].face_forces)
    end

    @testset "update_thermal_force! (deep crater changes the force; net force and torque re-summed)" begin
        path_obj = joinpath(@__DIR__, "shape", "icosahedron.obj")
        shape_plain = load_shape_obj(path_obj)
        shape_hier  = load_shape_obj(path_obj; as_hierarchical=true)
        add_roughness_models!(shape_hier, roughness_model)

        problem_plain = SingleAsteroidThermoPhysicalProblem(shape_plain, thermo_params, grid_params;
            with_self_shadowing=false, with_self_heating=false)
        problem_hier  = SingleAsteroidThermoPhysicalProblem(shape_hier, thermo_params, grid_params;
            with_self_shadowing=false, with_self_heating=false)
        state_plain = AsteroidThermoPhysicalModels._build_single_state(problem_plain, CrankNicolson())
        state_hier  = AsteroidThermoPhysicalModels._build_single_state(problem_hier,  CrankNicolson())
        AsteroidThermoPhysicalModels.init_temperature!(state_plain, 250.0)
        AsteroidThermoPhysicalModels.init_temperature!(state_hier,  250.0)

        # A few steps so that the crater walls develop a temperature contrast
        r☉ = SVector(0.3, -0.2, 0.9) * AsteroidThermoPhysicalModels.au2m
        for state in (state_plain, state_hier), _ in 1:5
            AsteroidThermoPhysicalModels.update_flux_sun!(state, r☉)
            AsteroidThermoPhysicalModels.update_flux_scat_single!(state)
            AsteroidThermoPhysicalModels.update_flux_rad_single!(state)
            AsteroidThermoPhysicalModels.update_temperature!(state, 100.0)
        end
        update_fluxes_and_force!(state_plain, r☉)
        update_fluxes_and_force!(state_hier,  r☉)

        @test all(state_hier.face_forces .!= state_plain.face_forces)
        @test all(f -> all(isfinite, f), state_hier.face_forces)
        @test state_hier.force  ≈ sum(state_hier.face_forces)
        @test state_hier.torque ≈ sum(r × f for (r, f) in zip(shape_hier.global_shape.face_centers, state_hier.face_forces))
    end

    @testset "update_thermal_force! (sky re-absorption term follows the global self-heating)" begin
        # On a globally concave shape, the photons that leave a roughness model towards the sky
        # can hit other global faces. Their momentum is added only when the global self-heating
        # is on; without it the face force is exactly the normalised sum of the sub-face forces.
        ẑ  = SVector(0.0, 0.0, 1.0)
        c₀ = AsteroidThermoPhysicalModels.c₀
        σ  = AsteroidThermoPhysicalModels.σ_SB
        R_vis, R_ir, ε = 0.1, 0.0, 0.9   # as in `thermo_params`
        r☉ = SVector(0.2, -0.1, 1.0) * AsteroidThermoPhysicalModels.au2m

        for with_self_heating in (true, false)
            shape_hier = create_shape_crater(0.4, 0.1; Nx=8, Ny=8, as_hierarchical=true)
            add_roughness_models!(shape_hier, create_shape_crater(0.4, 0.1; Nx=4, Ny=4))
            problem_hier = SingleAsteroidThermoPhysicalProblem(shape_hier, thermo_params, grid_params;
                with_self_shadowing=false, with_self_heating)
            state_hier = AsteroidThermoPhysicalModels._build_single_state(problem_hier, CrankNicolson())
            AsteroidThermoPhysicalModels.init_temperature!(state_hier, 250.0)
            update_fluxes_and_force!(state_hier, r☉)

            global_shape = shape_hier.global_shape
            n_sky = 0
            for (i, k) in enumerate(state_hier.face_roughness_indices)
                rs = state_hier.roughness_states[k]
                rshape = rs.problem.shape
                A_proj = sum(a * (n̂ ⋅ ẑ) for (a, n̂) in zip(rshape.face_areas, rshape.face_normals))
                patches_per_face = global_shape.face_areas[i] / A_proj
                patch = patches_per_face * transform_physical_vector_local_to_global(shape_hier, i, sum(rs.face_forces))

                if with_self_heating
                    P_sky = patches_per_face * sum(eachindex(rshape.faces)) do j
                        E_j   = R_vis * (rs.flux_sun[j] + rs.flux_scat[j]) + R_ir * rs.flux_rad[j] + ε * σ * rs.temperature[begin, j]^4
                        f_sky = max(0.0, 1 - sum(get_view_factors(rshape.face_visibility_graph, j)))
                        E_j * rshape.face_areas[j] * f_sky
                    end
                    sky = sum(zip(get_view_factors(global_shape.face_visibility_graph, i),
                                  get_visible_face_directions(global_shape.face_visibility_graph, i));
                              init=zero(SVector{3, Float64})) do (f, d̂)
                        P_sky / c₀ * f * d̂
                    end
                    norm(sky) > 0 && (n_sky += 1)
                    @test isapprox(state_hier.face_forces[i], patch + sky; rtol=1e-12)
                else
                    @test isapprox(state_hier.face_forces[i], patch; rtol=1e-12)
                end
            end
            with_self_heating && @test n_sky > 0
        end
    end
end
