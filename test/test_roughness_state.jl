#=
test_roughness_state.jl

Unit tests for SingleAsteroidThermoPhysicalState on a shape with surface roughness:
- _build_single_state for a ShapeModel carrying roughness models
- init_temperature! including the sub-face states (Real and AbstractMatrix)
- automatic preparation of the geometric data required for self-shadowing and self-heating
- update_flux_sun! on the global level, compared against the plain ShapeModel result
- update_flux_sun! on the sub-face level: gate on the global illumination, consistency with a
  standalone roughness model, and a near-flat roughness model reproducing the parent face
- update_flux_scat_single! / update_flux_rad_single! on the global level, likewise
- update_flux_scat_single! / update_flux_rad_single! on the sub-face level: self-heating inside
  the roughness model, the directional irradiation from the other global faces (against an
  independent reference implementation, on rough, flat and partly rough shapes), and its
  absence when the global self-heating is off
- update_temperature! on the global level, for all three solvers and for zero conductivity
- update_temperature! on the sub-face level: the sub-faces advance, a near-flat roughness model
  reproduces the temperature of its parent face, and zero conductivity gives radiative
  equilibrium on every sub-face
- update_thermal_force!: faces without roughness keep the plain result, a near-flat roughness
  model reproduces the force on its parent face, the roughness scale cancels, a deep crater
  changes the force and the net force/torque are re-summed, and the sky re-absorption term
  follows the global self-heating flag
=#

@testset "SingleAsteroidThermoPhysicalState with surface roughness" begin
    msg = """
    ┏━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━┓
    |         Test: TPM state with surface roughness         |
    ┗━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━┛
    """
    println(msg)

    # Build a ShapeModel from icosahedron with a crater roughness on all faces
    hier_shape = load_shape_obj(joinpath(@__DIR__, "shape", "icosahedron.obj"))
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

    n_global_faces = length(hier_shape.faces)

    @testset "problem construction" begin
        @test problem isa SingleAsteroidThermoPhysicalProblem
        @test problem.shape === hier_shape
        @test length(problem.thermo_params.conductivity) == n_global_faces
    end

    state = AsteroidThermoPhysicalModels._build_single_state(problem, CrankNicolson())

    @testset "_build_single_state returns a state with sub-face states" begin
        @test state isa AsteroidThermoPhysicalModels.SingleAsteroidThermoPhysicalState
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

    @testset "smooth shape has empty roughness vectors" begin
        shape_smooth = load_shape_obj(joinpath(@__DIR__, "shape", "icosahedron.obj"))
        problem_smooth = SingleAsteroidThermoPhysicalProblem(shape_smooth, thermo_params, grid_params;
            with_self_shadowing = false,
            with_self_heating   = false,
        )
        state_smooth = AsteroidThermoPhysicalModels._build_single_state(problem_smooth, CrankNicolson())

        @test isempty(state_smooth.face_roughness_indices)
        @test isempty(state_smooth.roughness_states)
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
        hier_shape2 = load_shape_obj(joinpath(@__DIR__, "shape", "icosahedron.obj"))
        add_roughness_models!(hier_shape2, roughness_model, 1)

        problem2 = SingleAsteroidThermoPhysicalProblem(hier_shape2, thermo_params, grid_params;
            with_self_shadowing = false,
            with_self_heating   = false,
        )
        state2 = AsteroidThermoPhysicalModels._build_single_state(problem2, CrankNicolson())

        n_faces2 = length(hier_shape2.faces)

        @test length(state2.face_roughness_indices) == n_faces2
        @test state2.face_roughness_indices[1] == 1       # face 1 has roughness
        @test all(state2.face_roughness_indices[2:end] .== 0)  # others don't
        @test length(state2.roughness_states) == 1
    end

    @testset "self-shadowing geometry prepared automatically" begin
        # Self-shadowing is evaluated on the (global) shape itself, so the visibility graph
        # and the maximum elevations must be built on it.
        hier_shape3 = load_shape_obj(joinpath(@__DIR__, "shape", "icosahedron.obj"))
        add_roughness_models!(hier_shape3, roughness_model)

        @test isnothing(hier_shape3.face_visibility_graph)
        @test isnothing(hier_shape3.face_max_elevations)

        SingleAsteroidThermoPhysicalProblem(hier_shape3, thermo_params, grid_params;
            with_self_shadowing = true,
            with_self_heating   = false,
        )

        @test !isnothing(hier_shape3.face_visibility_graph)
        @test !isnothing(hier_shape3.face_max_elevations)
    end

    @testset "self-heating geometry prepared automatically" begin
        hier_shape4 = load_shape_obj(joinpath(@__DIR__, "shape", "icosahedron.obj"))
        add_roughness_models!(hier_shape4, roughness_model)

        @test isnothing(hier_shape4.face_visibility_graph)

        SingleAsteroidThermoPhysicalProblem(hier_shape4, thermo_params, grid_params;
            with_self_shadowing = false,
            with_self_heating   = true,
        )

        @test !isnothing(hier_shape4.face_visibility_graph)
        @test isnothing(hier_shape4.face_max_elevations)  # not needed by self-heating
    end

    @testset "update_flux_sun! (global level)" begin
        # The global level of a state with roughness must reproduce the smooth-shape result
        # exactly: the sub-face machinery must not perturb the global solar flux.
        path_obj = joinpath(@__DIR__, "shape", "icosahedron.obj")
        au2m = AsteroidThermoPhysicalModels.au2m
        r☉s = [
            SVector(1.0, 0.0, 0.0) * au2m,        # 1.0 au
            SVector(0.3, -0.5, 0.8) * 1.2au2m,    # 1.2 au, oblique
        ]

        for with_self_shadowing in (false, true)
            shape_plain = load_shape_obj(path_obj)
            shape_hier  = load_shape_obj(path_obj)
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
                for (i, n̂) in enumerate(shape_hier.face_normals)
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

        shape_hier = load_shape_obj(path_obj)
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
        shape_hier = load_shape_obj(joinpath(@__DIR__, "shape", "icosahedron.obj"))
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
                shape_hier.face_normals[i] ⋅ r̂☉ < 1e-3 && continue  # grazing: skip
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
        shape_hier = load_shape_obj(joinpath(@__DIR__, "shape", "icosahedron.obj"))
        add_roughness_models!(shape_hier, create_shape_crater(0.4, 0.1; Nx=4, Ny=4), 1)
        problem_hier = SingleAsteroidThermoPhysicalProblem(shape_hier, thermo_params, grid_params;
            with_self_shadowing=false, with_self_heating=false)
        state_hier = AsteroidThermoPhysicalModels._build_single_state(problem_hier, CrankNicolson())

        # Point the Sun along the normal of face 1 so that it is certainly lit
        n̂₁ = shape_hier.face_normals[1]
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
        shape_hier = load_shape_obj(joinpath(@__DIR__, "shape", "icosahedron.obj"))
        add_roughness_models!(shape_hier, roughness_model)
        problem_hier = SingleAsteroidThermoPhysicalProblem(shape_hier, thermo_params, grid_params;
            with_self_shadowing=true, with_self_heating=false)
        state_hier = AsteroidThermoPhysicalModels._build_single_state(problem_hier, CrankNicolson())

        shape_hier.face_visibility_graph = nothing
        r☉ = SVector(1.0, 0.0, 0.0) * AsteroidThermoPhysicalModels.au2m
        @test_throws ArgumentError AsteroidThermoPhysicalModels.update_flux_sun!(state_hier, r☉)
    end

    @testset "update_flux_scat_single! / update_flux_rad_single! (global level)" begin
        # Self-heating vanishes identically on a convex shape, so the icosahedron used above
        # cannot tell a working implementation from a broken one. A crater is concave: its
        # faces see each other, and both self-heating terms are non-zero.
        make_crater() = create_shape_crater(0.4, 0.1; Nx=8, Ny=8)

        shape_plain = make_crater()
        shape_hier  = make_crater()
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

    # Reference implementation of the external irradiation of the sub-faces of global face `i`
    # (sub-state `k`), written independently of the package: for every visible global face
    # `j`, the emission of `j` towards `i` (directional for a rough `j`, Lambertian for a smooth
    # one) times the view factor, distributed over the sub-faces of `i` that see `j`. The
    # visibility is recomputed here with `update_illumination!`, so this also checks the
    # precomputed masks. `emission(state_or_substate, n)` gives the emitted quantity per face.
    function external_reference(state_hier, shape_hier, i, k, emission_global, emission_sub)
        graph = shape_hier.face_visibility_graph
        rs    = state_hier.roughness_states[k]
        model = rs.problem.shape
        ext   = zeros(length(model.faces))
        seen  = Vector{Bool}(undef, length(model.faces))
        total = 0.0     # Σ_j F_ij, the irradiance of the parent
        power = 0.0     # power received by the sub-faces, as from a distant source of each F_ij
        for (j, f, d̂) in zip(get_visible_face_indices(graph, i), get_view_factors(graph, i), get_visible_face_directions(graph, i))
            k_j = state_hier.face_roughness_indices[j]
            if k_j == 0
                E = emission_global(j)
            else
                rs_j = state_hier.roughness_states[k_j]; model_j = rs_j.problem.shape
                d̂_ji = transform_physical_vector_global_to_local(shape_hier, j, -d̂)
                seen_j = Vector{Bool}(undef, length(model_j.faces))
                update_illumination!(seen_j, model_j, d̂_ji; with_self_shadowing=true)
                E = sum(max(0.0, model_j.face_normals[n] ⋅ d̂_ji) * model_j.face_areas[n] * emission_sub(rs_j, n)
                        for n in eachindex(model_j.faces) if seen_j[n]; init=0.0) / (projected_area(model_j) * d̂_ji[3])
            end
            F = E * f
            total += F
            d̂_ij = transform_physical_vector_global_to_local(shape_hier, i, d̂)
            update_illumination!(seen, model, d̂_ij; with_self_shadowing=true)
            for m in eachindex(model.faces)
                seen[m] || continue
                c = model.face_normals[m] ⋅ d̂_ij
                c > 0 || continue
                ext[m] += F / d̂_ij[3] * c
                power  += F / d̂_ij[3] * c * model.face_areas[m]
            end
        end
        ext, total, power
    end
    ε_of(state) = state.problem.thermo_params.emissivity
    emission_rad_global(state)  = j -> ε_of(state)[j] * AsteroidThermoPhysicalModels.σ_SB * state.temperature[begin, j]^4
    emission_rad_sub            = (rs, n) -> ε_of(rs)[n] * AsteroidThermoPhysicalModels.σ_SB * rs.temperature[begin, n]^4
    emission_scat_global(state) = j -> state.problem.thermo_params.reflectance_vis[j] * state.flux_sun[j]
    emission_scat_sub           = (rs, n) -> rs.problem.thermo_params.reflectance_vis[n] * rs.flux_sun[n]

    @testset "update_flux_scat_single! / update_flux_rad_single! (sub-face level)" begin
        # Self-heating inside the roughness model must be exactly what a standalone state of
        # the roughness model gives, and the radiation the other global faces send towards the
        # sub-faces must be added on top — from the sub-faces of the neighbours' own roughness
        # models, direction by direction. A concave global crater makes the neighbours' terms
        # non-zero; a roughness crater makes the intra-model exchange non-zero.
        crater = create_shape_crater(0.4, 0.1; Nx=6, Ny=6)
        shape_hier = create_shape_crater(0.4, 0.1; Nx=8, Ny=8)
        add_roughness_models!(shape_hier, crater)
        problem_hier = SingleAsteroidThermoPhysicalProblem(shape_hier, thermo_params, grid_params;
            with_self_shadowing=false, with_self_heating=true)
        state_hier = AsteroidThermoPhysicalModels._build_single_state(problem_hier, CrankNicolson())

        # Self-heating inside the roughness model is always on
        @test all(rs.problem.with_self_heating for rs in state_hier.roughness_states)
        @test !isnothing(crater.face_visibility_graph)

        # The pair visibility exists for every roughness state, with one mask per visible face
        @test length(state_hier.roughness_neighbours) == length(state_hier.roughness_states)
        for (i, k) in enumerate(state_hier.face_roughness_indices)
            nb = state_hier.roughness_neighbours[k]
            @test length(nb.visible_sub_faces) == n_visible_faces(shape_hier.face_visibility_graph, i)
            @test all(length(v) == length(crater.faces) for v in nb.visible_sub_faces)
            # Every neighbour carries roughness here, so every reverse position points back at i
            for (p, j) in enumerate(get_visible_face_indices(shape_hier.face_visibility_graph, i))
                @test get_visible_face_indices(shape_hier.face_visibility_graph, j)[nb.position_in_neighbour[p]] == i
            end
        end

        # Standalone state of the same roughness model, with the same flags as the sub-states
        problem_alone = SingleAsteroidThermoPhysicalProblem(crater, thermo_params, grid_params;
            with_self_shadowing=true, with_self_heating=true)
        state_alone = AsteroidThermoPhysicalModels._build_single_state(problem_alone, CrankNicolson())

        # Vary the temperature from face to face so the radiation term is not degenerate
        n_faces = length(shape_hier.faces)
        T₀ = repeat(reshape(range(200.0, 300.0; length=n_faces) |> collect, 1, n_faces), grid_params.n_depth)
        AsteroidThermoPhysicalModels.init_temperature!(state_hier, T₀)

        r☉ = SVector(0.2, -0.1, 1.0) * AsteroidThermoPhysicalModels.au2m
        AsteroidThermoPhysicalModels.update_flux_sun!(state_hier, r☉)
        AsteroidThermoPhysicalModels.update_flux_scat_single!(state_hier)
        AsteroidThermoPhysicalModels.update_flux_rad_single!(state_hier)

        # The global fluxes are non-zero somewhere, so the shape really is concave
        @test any(>(0), state_hier.flux_scat)
        @test any(>(0), state_hier.flux_rad)

        n_checked = n_irradiated = 0
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

            ext_rad,  total_rad,  power_rad  = external_reference(state_hier, shape_hier, i, k, emission_rad_global(state_hier),  emission_rad_sub)
            ext_scat, total_scat, power_scat = external_reference(state_hier, shape_hier, i, k, emission_scat_global(state_hier), emission_scat_sub)
            @test rs.flux_rad  ≈ state_alone.flux_rad  .+ ext_rad
            @test rs.flux_scat ≈ state_alone.flux_scat .+ ext_scat

            # Each neighbour irradiates the sub-faces like a distant source of irradiance F_ij
            # (the same rule as the sunlight): the power received is F_ij / cos θ times the
            # projected area of the sub-faces seen from that direction. This equals F_ij A_proj
            # up to the shadowing discretisation of the sub-faces (a sub-face is seen or not
            # from the ray through its centre), which is a few percent on this coarse crater.
            model = rs.problem.shape
            @test sum((rs.flux_rad .- state_alone.flux_rad) .* model.face_areas) ≈ power_rad
            total_rad > 0 && (n_irradiated += 1)   # the flat rim of the global crater sees nothing

            # ...and it is directional: the sub-faces that see no neighbour receive nothing,
            # the others do, unlike an isotropic (sky-view-factor) distribution
            seen_any = falses(length(model.faces))
            for d̂ in get_visible_face_directions(shape_hier.face_visibility_graph, i)
                d̂_ij = transform_physical_vector_global_to_local(shape_hier, i, d̂)
                seen = Vector{Bool}(undef, length(model.faces))
                update_illumination!(seen, model, d̂_ij; with_self_shadowing=true)
                seen_any .|= seen .& (map(n̂ -> n̂ ⋅ d̂_ij, model.face_normals) .> 0)
            end
            @test all(iszero, ext_rad[.!seen_any])
            @test all(>(0), ext_rad[seen_any])
        end
        @test n_checked > 0
        @test n_irradiated > 0
    end

    @testset "external irradiation differs from an isotropic distribution on a deep crater" begin
        # With a deep roughness crater, the walls facing a neighbour receive its radiation and the
        # walls behind do not; distributing the same power by the sky view factor (the isotropic
        # limit) gives a different pattern. On a flat model the two coincide (next testset).
        crater = create_shape_crater(0.5, 0.3; Nx=6, Ny=6)
        shape_hier = create_shape_crater(0.4, 0.1; Nx=8, Ny=8)
        add_roughness_models!(shape_hier, crater)
        problem_hier = SingleAsteroidThermoPhysicalProblem(shape_hier, thermo_params, grid_params;
            with_self_shadowing=false, with_self_heating=true)
        state_hier = AsteroidThermoPhysicalModels._build_single_state(problem_hier, CrankNicolson())
        AsteroidThermoPhysicalModels.init_temperature!(state_hier, 250.0)

        r☉ = SVector(0.2, -0.1, 1.0) * AsteroidThermoPhysicalModels.au2m
        AsteroidThermoPhysicalModels.update_flux_sun!(state_hier, r☉)
        AsteroidThermoPhysicalModels.update_flux_scat_single!(state_hier)
        AsteroidThermoPhysicalModels.update_flux_rad_single!(state_hier)

        f_sky = [max(0.0, 1 - sum(get_view_factors(crater.face_visibility_graph, m))) for m in eachindex(crater.faces)]
        n_differ = 0
        for (i, k) in enumerate(state_hier.face_roughness_indices)
            ext_rad, total_rad, _ = external_reference(state_hier, shape_hier, i, k, emission_rad_global(state_hier), emission_rad_sub)
            total_rad > 0 || continue
            isotropic = f_sky .* total_rad     # same power, spread as the sky-view-factor rule would
            isapprox(ext_rad, isotropic; rtol=1e-3) || (n_differ += 1)
        end
        @test n_differ > 0
    end

    @testset "global self-heating off removes the external term only" begin
        # With `with_self_heating = false` on the global problem, no radiation is exchanged
        # between global faces, no pair visibility is built, and the sub-faces receive nothing
        # from outside — but the self-heating inside each roughness model is still computed.
        # This is what makes an on/off comparison of the global self-heating isolate the
        # external contribution.
        crater = create_shape_crater(0.4, 0.1; Nx=6, Ny=6)
        shape_hier = create_shape_crater(0.4, 0.1; Nx=8, Ny=8)
        add_roughness_models!(shape_hier, crater)
        problem_hier = SingleAsteroidThermoPhysicalProblem(shape_hier, thermo_params, grid_params;
            with_self_shadowing=false, with_self_heating=false)
        state_hier = AsteroidThermoPhysicalModels._build_single_state(problem_hier, CrankNicolson())
        AsteroidThermoPhysicalModels.init_temperature!(state_hier, 250.0)
        @test isempty(state_hier.roughness_neighbours)

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
        # A flat roughness model has no walls: every sub-face sees every neighbour with the
        # inclination of the parent, and the neighbours (flat too) radiate as their smooth
        # faces. The directional exchange therefore reduces to the smooth-surface flux of the
        # parent face — the isotropic limit — and each sub-face must end up with exactly it.
        shape_hier = create_shape_crater(0.4, 0.1; Nx=8, Ny=8)
        add_roughness_models!(shape_hier, create_shape_crater(0.4, 1e-6; Nx=4, Ny=4))
        problem_hier = SingleAsteroidThermoPhysicalProblem(shape_hier, thermo_params, grid_params;
            with_self_shadowing=false, with_self_heating=true)
        state_hier = AsteroidThermoPhysicalModels._build_single_state(problem_hier, CrankNicolson())
        AsteroidThermoPhysicalModels.init_temperature!(state_hier, 250.0)

        r☉ = SVector(0.2, -0.1, 1.0) * AsteroidThermoPhysicalModels.au2m
        AsteroidThermoPhysicalModels.update_flux_sun!(state_hier, r☉)
        AsteroidThermoPhysicalModels.update_flux_scat_single!(state_hier)
        AsteroidThermoPhysicalModels.update_flux_rad_single!(state_hier)

        # A model flat to 1e-6 is not exactly flat: the tilt of its normals (~1e-5) is felt
        # for the neighbours seen at grazing angles (cos θ down to ~1e-3 on this crater), and
        # the solar flux of its sub-faces differs from the parent's by ~1e-5 of the solar
        # constant, which enters the reflected sunlight of the neighbours. Compare with an
        # absolute tolerance scaled by the largest parent flux: a relative one would fail on
        # parents whose flux is exactly zero.
        atol_scat = 1e-4 * maximum(state_hier.flux_scat)
        atol_rad  = 1e-4 * maximum(state_hier.flux_rad)
        @test atol_scat > 0 && atol_rad > 0
        for (i, k) in enumerate(state_hier.face_roughness_indices)
            rs = state_hier.roughness_states[k]
            @test all(isapprox.(rs.flux_scat, state_hier.flux_scat[i]; atol=atol_scat))
            @test all(isapprox.(rs.flux_rad,  state_hier.flux_rad[i];  atol=atol_rad))
        end
    end

    @testset "external irradiation with roughness on part of the faces" begin
        # Neighbours without a roughness model radiate as smooth Lambertian faces; the rough
        # faces among them radiate directionally. The reference implementation handles both.
        crater = create_shape_crater(0.4, 0.1; Nx=6, Ny=6)
        shape_hier = create_shape_crater(0.4, 0.1; Nx=8, Ny=8)
        rough_faces = 1:2:length(shape_hier.faces)
        for i in rough_faces
            add_roughness_models!(shape_hier, crater, i)
        end
        problem_hier = SingleAsteroidThermoPhysicalProblem(shape_hier, thermo_params, grid_params;
            with_self_shadowing=false, with_self_heating=true)
        state_hier = AsteroidThermoPhysicalModels._build_single_state(problem_hier, CrankNicolson())
        n_faces = length(shape_hier.faces)
        T₀ = repeat(reshape(range(200.0, 300.0; length=n_faces) |> collect, 1, n_faces), grid_params.n_depth)
        AsteroidThermoPhysicalModels.init_temperature!(state_hier, T₀)

        r☉ = SVector(0.2, -0.1, 1.0) * AsteroidThermoPhysicalModels.au2m
        AsteroidThermoPhysicalModels.update_flux_sun!(state_hier, r☉)
        AsteroidThermoPhysicalModels.update_flux_scat_single!(state_hier)
        AsteroidThermoPhysicalModels.update_flux_rad_single!(state_hier)

        problem_alone = SingleAsteroidThermoPhysicalProblem(crater, thermo_params, grid_params;
            with_self_shadowing=true, with_self_heating=true)
        state_alone = AsteroidThermoPhysicalModels._build_single_state(problem_alone, CrankNicolson())

        n_mixed = 0
        for (i, k) in enumerate(state_hier.face_roughness_indices)
            k == 0 && continue
            rs = state_hier.roughness_states[k]
            nb = state_hier.roughness_neighbours[k]
            js = get_visible_face_indices(shape_hier.face_visibility_graph, i)
            any(j -> state_hier.face_roughness_indices[j] == 0, js) && any(j -> state_hier.face_roughness_indices[j] != 0, js) && (n_mixed += 1)
            # Reverse positions: 0 for smooth neighbours, the position of i otherwise
            for (p, j) in enumerate(js)
                if state_hier.face_roughness_indices[j] == 0
                    @test nb.position_in_neighbour[p] == 0
                else
                    @test get_visible_face_indices(shape_hier.face_visibility_graph, j)[nb.position_in_neighbour[p]] == i
                end
            end

            AsteroidThermoPhysicalModels.init_temperature!(state_alone, T₀[begin, i])
            AsteroidThermoPhysicalModels.update_flux_sun!(state_alone, transform_physical_vector_global_to_local(shape_hier, i, r☉))
            AsteroidThermoPhysicalModels.update_flux_scat_single!(state_alone)
            AsteroidThermoPhysicalModels.update_flux_rad_single!(state_alone)
            ext_rad, _, _  = external_reference(state_hier, shape_hier, i, k, emission_rad_global(state_hier),  emission_rad_sub)
            ext_scat, _, _ = external_reference(state_hier, shape_hier, i, k, emission_scat_global(state_hier), emission_scat_sub)
            @test rs.flux_rad  ≈ state_alone.flux_rad  .+ ext_rad
            @test rs.flux_scat ≈ state_alone.flux_scat .+ ext_scat
            @test all(isfinite, rs.flux_rad)
        end
        @test n_mixed > 0
    end

    @testset "self-heating disabled leaves the global fluxes at zero" begin
        shape_hier = create_shape_crater(0.4, 0.1; Nx=8, Ny=8)
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
        # The global faces of a state with roughness are solved independently of their roughness
        # models, so after any number of steps they must match the plain ShapeModel exactly.
        # A concave crater with self-heating exercises every flux term in the surface balance.
        make_crater() = create_shape_crater(0.4, 0.1; Nx=8, Ny=8)

        n_faces = length(make_crater().faces)
        T₀ = repeat(reshape(range(200.0, 300.0; length=n_faces) |> collect, 1, n_faces),
                    grid_params.n_depth)
        r☉ = SVector(0.2, -0.1, 1.0) * AsteroidThermoPhysicalModels.au2m
        Δt = 100.0  # λ = αΔt/Δz² ≈ 0.12 keeps the explicit Euler step stable
        n_steps = 10

        for algorithm in (CrankNicolson(), ImplicitEuler(), ExplicitEuler())
            shape_plain = make_crater()
            shape_hier  = make_crater()
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
            shape_hier = load_shape_obj(path_obj)
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
                state_hier.illuminated_faces[i] && shape_hier.face_normals[i] ⋅ r̂☉ < 0.05 && continue
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
        shape_hier  = create_shape_crater(0.4, 0.1; Nx=8, Ny=8)
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
        shape_hier  = load_shape_obj(path_obj)
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
        shape_hier  = load_shape_obj(path_obj)
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
            shape_hier = load_shape_obj(path_obj)
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
        shape_hier  = load_shape_obj(path_obj)
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
        @test state_hier.torque ≈ sum(r × f for (r, f) in zip(shape_hier.face_centers, state_hier.face_forces))
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
            shape_hier = create_shape_crater(0.4, 0.1; Nx=8, Ny=8)
            add_roughness_models!(shape_hier, create_shape_crater(0.4, 0.1; Nx=4, Ny=4))
            problem_hier = SingleAsteroidThermoPhysicalProblem(shape_hier, thermo_params, grid_params;
                with_self_shadowing=false, with_self_heating)
            state_hier = AsteroidThermoPhysicalModels._build_single_state(problem_hier, CrankNicolson())
            AsteroidThermoPhysicalModels.init_temperature!(state_hier, 250.0)
            update_fluxes_and_force!(state_hier, r☉)

            n_sky = 0
            for (i, k) in enumerate(state_hier.face_roughness_indices)
                rs = state_hier.roughness_states[k]
                rshape = rs.problem.shape
                A_proj = sum(a * (n̂ ⋅ ẑ) for (a, n̂) in zip(rshape.face_areas, rshape.face_normals))
                patches_per_face = shape_hier.face_areas[i] / A_proj
                patch = patches_per_face * transform_physical_vector_local_to_global(shape_hier, i, sum(rs.face_forces))

                if with_self_heating
                    P_sky = patches_per_face * sum(eachindex(rshape.faces)) do j
                        E_j   = R_vis * (rs.flux_sun[j] + rs.flux_scat[j]) + R_ir * rs.flux_rad[j] + ε * σ * rs.temperature[begin, j]^4
                        f_sky = max(0.0, 1 - sum(get_view_factors(rshape.face_visibility_graph, j)))
                        E_j * rshape.face_areas[j] * f_sky
                    end
                    sky = sum(zip(get_view_factors(shape_hier.face_visibility_graph, i),
                                  get_visible_face_directions(shape_hier.face_visibility_graph, i));
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
