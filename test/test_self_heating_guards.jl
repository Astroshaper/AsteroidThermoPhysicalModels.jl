#=
test_self_heating_guards.jl

The terms that depend on the face visibility graph must fail loudly when the graph is missing
and the flag that needs it is enabled, and the re-absorption recoil term must follow
`with_self_heating` rather than the mere presence of the graph.
=#

@testset "Self-heating guards" begin
    msg = """
    ┏━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━┓
    |              Test: Self-heating guards                 |
    ┗━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━┛
    """
    println(msg)

    thermo_params = ThermoParams(
        conductivity    = 0.1,
        density         = 1000.0,
        heat_capacity   = 700.0,
        reflectance_vis = 0.1,
        reflectance_ir  = 0.0,
        emissivity      = 0.9,
    )
    grid_params = GridParams(; z_max=0.1, n_depth=10)
    r☉ = SVector(0.2, -0.1, 1.0) * AsteroidThermoPhysicalModels.au2m

    # A concave crater: its faces see each other, so the terms under test are non-zero
    function make_state(; with_self_shadowing, with_self_heating)
        shape   = create_shape_crater(0.4, 0.1; Nx=8, Ny=8)
        problem = SingleAsteroidThermoPhysicalProblem(shape, thermo_params, grid_params;
            with_self_shadowing, with_self_heating)
        state = AsteroidThermoPhysicalModels._build_single_state(problem, CrankNicolson())
        AsteroidThermoPhysicalModels.init_temperature!(state, 250.0)
        state
    end

    # Direct recoil only: -2/3 E A / c n̂, with E from the fluxes the state currently holds
    function direct_recoil(state, i)
        shape = state.problem.shape
        R_vis = state.problem.thermo_params.reflectance_vis[i]
        R_ir  = state.problem.thermo_params.reflectance_ir[i]
        ε     = state.problem.thermo_params.emissivity[i]
        Eᵢ = R_vis * state.flux_sun[i] + R_vis * state.flux_scat[i] + R_ir * state.flux_rad[i] +
             ε * AsteroidThermoPhysicalModels.σ_SB * state.temperature[begin, i]^4
        -2/3 * Eᵢ * shape.face_areas[i] / AsteroidThermoPhysicalModels.c₀ * shape.face_normals[i]
    end

    @testset "scattering and radiation raise without the visibility graph" begin
        state = make_state(with_self_shadowing=false, with_self_heating=true)
        AsteroidThermoPhysicalModels.update_flux_sun!(state, r☉)
        state.problem.shape.face_visibility_graph = nothing
        @test_throws ErrorException AsteroidThermoPhysicalModels.update_flux_scat_single!(state)
        @test_throws ErrorException AsteroidThermoPhysicalModels.update_flux_rad_single!(state)
    end

    @testset "scattering and radiation stay silent when self-heating is off" begin
        # No graph and no flag: nothing to compute, nothing to complain about
        state = make_state(with_self_shadowing=false, with_self_heating=false)
        @test isnothing(state.problem.shape.face_visibility_graph)
        AsteroidThermoPhysicalModels.update_flux_sun!(state, r☉)
        AsteroidThermoPhysicalModels.update_flux_scat_single!(state)
        AsteroidThermoPhysicalModels.update_flux_rad_single!(state)
        @test all(state.flux_scat .== 0.0)
        @test all(state.flux_rad  .== 0.0)
    end

    @testset "re-absorption recoil follows with_self_heating, not the graph" begin
        # Self-shadowing builds the graph while self-heating is off: the graph is present, but
        # the re-absorption term must not be applied, so each face force is the direct recoil.
        state = make_state(with_self_shadowing=true, with_self_heating=false)
        @test !isnothing(state.problem.shape.face_visibility_graph)
        AsteroidThermoPhysicalModels.update_flux_all!(state, r☉)
        AsteroidThermoPhysicalModels.update_thermal_force!(state)
        for i in eachindex(state.problem.shape.faces)
            @test state.face_forces[i] ≈ direct_recoil(state, i)
        end

        # With self-heating on, the re-absorption term is applied: on a concave shape at least
        # one face force departs from the direct recoil
        state_sh = make_state(with_self_shadowing=true, with_self_heating=true)
        AsteroidThermoPhysicalModels.update_flux_all!(state_sh, r☉)
        AsteroidThermoPhysicalModels.update_thermal_force!(state_sh)
        @test any(i -> !(state_sh.face_forces[i] ≈ direct_recoil(state_sh, i)),
                  eachindex(state_sh.problem.shape.faces))
    end

    @testset "thermal force raises without the visibility graph when self-heating is on" begin
        state = make_state(with_self_shadowing=false, with_self_heating=true)
        AsteroidThermoPhysicalModels.update_flux_all!(state, r☉)
        state.problem.shape.face_visibility_graph = nothing
        @test_throws ErrorException AsteroidThermoPhysicalModels.update_thermal_force!(state)
    end
end
