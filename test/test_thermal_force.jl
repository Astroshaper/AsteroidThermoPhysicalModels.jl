#=
test_thermal_force.jl

The net thermal force is the plain sum of the face forces and the net torque the sum of
their moments about the origin. Momentum conservation pins the first: on an isothermal
closed body every photon leaving is balanced by one leaving elsewhere, so the net recoil
must vanish exactly.
=#

@testset "Thermal force" begin
    msg = """
    ┏━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━┓
    |                  Test: Thermal force                   |
    ┗━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━┛
    """
    println(msg)

    # A closed, irregular body: face normals are not radial, so the sum of the face forces
    # and the sum of their radial projections differ. Convex enough that no visibility
    # graph is needed for what is tested here.
    shape = load_shape_obj(joinpath(@__DIR__, "shape", "ryugu_test.obj"); scale=1000)
    thermo_params = ThermoParams(
        conductivity    = 0.1,
        density         = 1000.0,
        heat_capacity   = 700.0,
        reflectance_vis = 0.1,
        reflectance_ir  = 0.0,
        emissivity      = 0.9,
    )
    grid_params = GridParams(; z_max=0.1, n_depth=10)
    problem = SingleAsteroidThermoPhysicalProblem(shape, thermo_params, grid_params;
        with_self_shadowing=false, with_self_heating=false)
    state = AsteroidThermoPhysicalModels._build_single_state(problem, CrankNicolson())
    c₀ = AsteroidThermoPhysicalModels.c₀

    @testset "isothermal closed body has no net recoil" begin
        # No sunlight, uniform temperature: every face emits the same E, so Fᵢ ∝ −aᵢ n̂ᵢ and
        # Σᵢ aᵢ n̂ᵢ = 0 on a closed surface. The torque of a uniform normal pressure about
        # any point vanishes as well.
        T = 250.0
        AsteroidThermoPhysicalModels.init_temperature!(state, T)
        AsteroidThermoPhysicalModels.update_thermal_force!(state)

        E = 0.9 * AsteroidThermoPhysicalModels.σ_SB * T^4
        F_scale = 2/3 * E * sum(shape.face_areas) / c₀          # one-sided force scale
        R = maximum(norm, shape.face_centers)

        @test norm(state.force)  < 1e-12 * F_scale
        @test norm(state.torque) < 1e-12 * F_scale * R

        # The radial projection that used to be summed does not vanish on this body — this
        # is what the test guards against
        F_proj = sum(normalize(shape.face_centers[i]) ⋅ state.face_forces[i] * normalize(shape.face_centers[i])
                     for i in eachindex(shape.faces))
        @test norm(F_proj) > 1e-4 * F_scale
    end

    @testset "net force and torque are the sums over the faces" begin
        # A temperature gradient makes the recoil asymmetric and the net force non-zero
        n_faces = length(shape.faces)
        T₀ = repeat(reshape(range(200.0, 300.0; length=n_faces) |> collect, 1, n_faces), grid_params.n_depth)
        AsteroidThermoPhysicalModels.init_temperature!(state, T₀)
        AsteroidThermoPhysicalModels.update_thermal_force!(state)

        F_sum = sum(state.face_forces)
        τ_sum = sum(shape.face_centers[i] × state.face_forces[i] for i in eachindex(shape.faces))
        @test norm(F_sum) > 0
        @test state.force  ≈ F_sum
        @test state.torque ≈ τ_sum
    end
end
