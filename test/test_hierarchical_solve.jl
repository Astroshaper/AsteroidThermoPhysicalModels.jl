#=
test_hierarchical_solve.jl

End-to-end tests of `solve` on a `HierarchicalShapeModel`:
- with roughness on one face only, every global-level output equals the plain ShapeModel run
  exactly, and the force and power of that face are counted from its crater
- a near-flat roughness model on every face reproduces the plain net force, torque and power
- the energy balance closes with roughness faces counted as representative patches
- export_solution writes the usual files, and show_progress runs
=#

@testset "solve on HierarchicalShapeModel" begin
    msg = """
    ┏━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━┓
    |          Test: solve on HierarchicalShapeModel         |
    ┗━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━┛
    """
    println(msg)

    au2m = AsteroidThermoPhysicalModels.au2m
    σ_SB = AsteroidThermoPhysicalModels.σ_SB
    c₀   = AsteroidThermoPhysicalModels.c₀

    P            = 8 * 3600.0  # Rotation period [s]
    n_step_cycle = 36
    path_obj     = joinpath(@__DIR__, "shape", "icosahedron.obj")

    # Rotation about the z axis by `θ`
    rotz(θ) = SMatrix{3,3,Float64}(cos(θ), sin(θ), 0, -sin(θ), cos(θ), 0, 0, 0, 1)

    thermo_params = ThermoParams(
        conductivity    = 0.1,
        density         = 1000.0,
        heat_capacity   = 600.0,
        reflectance_vis = 0.1,
        reflectance_ir  = 0.0,
        emissivity      = 0.9,
    )
    grid_params = GridParams(; z_max=0.1, n_depth=5)

    roughness_model = create_shape_crater(0.4, 0.1; Nx=4, Ny=4)

    # `n_cycle` rotations about the body z axis with the Sun fixed on the inertial x axis at 1 au
    function make_ephem(n_cycle; with_rotation)
        et_range = range(0.0, P * n_cycle; length=n_step_cycle * n_cycle + 1)
        times    = collect(et_range)
        r_sun    = [rotz(-2π * et / P) * SVector(au2m, 0.0, 0.0) for et in et_range]
        if with_rotation
            R_b2i = [rotz(2π * et / P) for et in et_range]
            return times, SingleAsteroidEphemerides(times, r_sun, R_b2i)
        else
            return times, SingleAsteroidEphemerides(times, r_sun)
        end
    end

    make_problem(shape; with_self_heating=false) =
        SingleAsteroidThermoPhysicalProblem(shape, thermo_params, grid_params;
            with_self_shadowing=false, with_self_heating)

    @testset "partial roughness: global-level outputs equal the plain run" begin
        times, ephem = make_ephem(2; with_rotation=false)
        output_times = times[end-n_step_cycle:end]
        output = SingleAsteroidOutputSpec(output_times, [1, 7];
            save_surface_temperature    = true,
            save_subsurface_temperature = true,
            save_face_forces            = true,
        )

        shape_plain = load_shape_obj(path_obj)
        shape_hier  = load_shape_obj(path_obj; as_hierarchical=true)
        add_roughness_models!(shape_hier, roughness_model, 1)

        sol_plain = solve(make_problem(shape_plain), CrankNicolson(); ephem, output, initial_temperature=200.0)
        sol_hier  = solve(make_problem(shape_hier),  CrankNicolson(); ephem, output, initial_temperature=200.0)

        @test sol_hier isa SingleAsteroidThermoPhysicalSolution

        # The global faces are solved independently of their roughness models, so every
        # global-level temperature is bit-identical to the plain run — face 1 included.
        @test sol_hier.surface_temperature       == sol_plain.surface_temperature
        @test sol_hier.subsurface_temperature[1] == sol_plain.subsurface_temperature[1]
        @test sol_hier.subsurface_temperature[7] == sol_plain.subsurface_temperature[7]

        # Faces without a roughness model keep the plain force; face 1 carries the crater's.
        @test sol_hier.face_forces[2:end, :] == sol_plain.face_forces[2:end, :]
        @test all(sol_hier.face_forces[1, :] .!= sol_plain.face_forces[1, :])

        # Face 1 is counted from its sub-faces in the power budget.
        @test sol_hier.absorbed_power != sol_plain.absorbed_power
        @test sol_hier.emitted_power  != sol_plain.emitted_power
        @test all(isfinite, sol_hier.absorbed_power)
        @test all(isfinite, sol_hier.emitted_power)
    end

    @testset "near-flat roughness reproduces the plain net force, torque and power" begin
        times, ephem = make_ephem(2; with_rotation=true)
        output_times = times[end-n_step_cycle:end]
        output = SingleAsteroidOutputSpec(output_times, Int[];
            save_surface_temperature    = false,
            save_subsurface_temperature = false,
            save_face_forces            = false,
            save_forces                 = true,
            save_torques                = true,
        )

        shape_plain = load_shape_obj(path_obj)
        shape_hier  = load_shape_obj(path_obj; as_hierarchical=true)
        add_roughness_models!(shape_hier, create_shape_crater(0.4, 1e-6; Nx=4, Ny=4))

        sol_plain = solve(make_problem(shape_plain), CrankNicolson(); ephem, output, initial_temperature=200.0)
        sol_hier  = solve(make_problem(shape_hier),  CrankNicolson(); ephem, output, initial_temperature=200.0)

        # The net force nearly cancels on a symmetric body, so compare against the scale of the
        # one-sided recoil rather than relatively: the 1e-6 tilt of the sub-face normals is an
        # absolute error on each face force.
        A_total = sum(shape_plain.face_areas)
        r_max   = maximum(norm, shape_plain.face_centers)
        F_scale = 2/3 * 0.9 * σ_SB * 300.0^4 * A_total / c₀
        @test all(isapprox.(sol_hier.forces,  sol_plain.forces;  atol=1e-4 * F_scale))
        @test all(isapprox.(sol_hier.torques, sol_plain.torques; atol=1e-4 * F_scale * r_max))
        @test any(F -> norm(F) > 1e-3 * F_scale, sol_plain.forces)  # the comparison is not vacuous

        P_scale = maximum(sol_plain.absorbed_power)
        @test all(isapprox.(sol_hier.absorbed_power, sol_plain.absorbed_power; atol=1e-4 * P_scale))
        @test all(isapprox.(sol_hier.emitted_power,  sol_plain.emitted_power;  atol=1e-4 * P_scale))
    end

    @testset "energy balance closes with roughness faces as representative patches" begin
        n_cycle = 10
        times, ephem = make_ephem(n_cycle; with_rotation=false)
        output = SingleAsteroidOutputSpec(times[end:end], Int[];
            save_surface_temperature    = true,
            save_subsurface_temperature = false,
        )

        shape_hier = load_shape_obj(path_obj; as_hierarchical=true)
        add_roughness_models!(shape_hier, roughness_model)

        sol = solve(make_problem(shape_hier; with_self_heating=true), CrankNicolson();
            ephem, output, initial_temperature=200.0)

        last_cycle = length(times)-n_step_cycle+1:length(times)
        ratio = sum(sol.emitted_power[last_cycle]) / sum(sol.absorbed_power[last_cycle])
        @test abs(ratio - 1) < 0.05
    end

    @testset "export_solution and show_progress" begin
        times, ephem = make_ephem(1; with_rotation=true)
        output_times = times[end-5:end]
        output = SingleAsteroidOutputSpec(output_times, [1];
            save_surface_temperature    = true,
            save_subsurface_temperature = true,
            save_face_forces            = true,
            save_forces                 = true,
            save_torques                = true,
        )

        shape_hier = load_shape_obj(path_obj; as_hierarchical=true)
        add_roughness_models!(shape_hier, roughness_model)

        sol = solve(make_problem(shape_hier), CrankNicolson();
            ephem, output, initial_temperature=200.0, show_progress=true)

        dir = mktempdir()
        export_solution(dir, sol)
        for file in ("diagnostics.csv", "surface_temperature.csv", "subsurface_temperature.csv",
                     "thermal_face_forces.csv", "thermal_net_forces.csv")
            @test isfile(joinpath(dir, file))
        end
        # One header line plus one row per output time
        @test countlines(joinpath(dir, "thermal_net_forces.csv")) == length(output_times) + 1
    end
end
