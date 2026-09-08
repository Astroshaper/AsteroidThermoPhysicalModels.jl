#=
test_directional_radiance.jl

Direction-dependent radiance and brightness temperature of rough facets:
- a near-flat roughness model reproduces the smooth Lambertian facet in every direction
- an isothermal crater is Lambertian (normalisation by the projected area)
- a sunlit crater beams: brighter towards the Sun than away from it
- the Planck inversion round-trips
=#

@testset "Directional radiance" begin
    msg = """
    ┏━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━┓
    |              Test: Directional radiance                |
    ┗━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━┛
    """
    println(msg)

    ATPM = AsteroidThermoPhysicalModels
    au2m = ATPM.au2m
    path_obj = joinpath(@__DIR__, "shape", "icosahedron.obj")
    ε = 0.9
    thermo_params = ThermoParams(conductivity=0.1, density=1000.0, heat_capacity=600.0,
                                 reflectance_vis=0.1, reflectance_ir=0.0, emissivity=ε)
    grid_params = GridParams(; z_max=0.1, n_depth=5)
    λ_tir = 10e-6

    # Sun fixed on the body +x axis at 1 au, no rotation: the temperature field settles
    # towards a steady day/night pattern
    P = 8 * 3600.0
    times = collect(range(0.0, P; length=73))
    ephem = SingleAsteroidEphemerides(times, [SVector(au2m, 0.0, 0.0) for _ in times])
    output_times = times[end:end]

    function run(roughness_model)
        shape = load_shape_obj(path_obj; as_hierarchical=true)
        add_roughness_models!(shape, roughness_model)
        n_face = length(shape.global_shape.faces)
        problem = SingleAsteroidThermoPhysicalProblem(shape, thermo_params, grid_params;
            with_self_shadowing=false, with_self_heating=false)
        output = SingleAsteroidOutputSpec(output_times, Int[];
            save_subsurface_temperature=false,
            roughness_face_ids=collect(1:n_face), save_roughness_surface_temperature=true)
        sol = solve(problem, CrankNicolson(); ephem, output, initial_temperature=200.0, show_progress=false)
        shape, problem, sol
    end

    @testset "near-flat roughness is Lambertian in every direction" begin
        shape, problem, sol = run(create_shape_crater(0.4, 1e-6; Nx=4, Ny=4))
        T = sol.surface_temperature[:, 1]
        for d̂ in (SVector(1.0, 0.0, 0.0), SVector(0.3, -0.5, 0.8), SVector(-0.2, 0.9, 0.1))
            d̂ = normalize(d̂)
            L   = directional_radiance(problem, sol, 1, d̂)
            T_b = brightness_temperature(problem, sol, 1, d̂)
            L_λ = directional_radiance(problem, sol, 1, d̂; λ=λ_tir)
            n_seen = 0
            for i in eachindex(T)
                cosθ = shape.global_shape.face_normals[i] ⋅ d̂
                if cosθ > 0.05
                    n_seen += 1
                    @test L[i]   ≈ ε * ATPM.σ_SB * T[i]^4 / π                      rtol=1e-4
                    @test T_b[i] ≈ ε^(1/4) * T[i]                                   rtol=1e-4
                    @test L_λ[i] ≈ ε * ATPM.blackbody_radiance(λ_tir, T[i]) / π      rtol=1e-4
                elseif cosθ < 0
                    @test isnan(L[i]) && isnan(T_b[i]) && isnan(L_λ[i])
                end
            end
            @test n_seen > 0
        end
    end

    @testset "isothermal crater is Lambertian" begin
        # Every sub-facet at the same temperature: the visible emission per projected area is
        # ε B(T)/π regardless of direction. Exact at normal incidence; at oblique angles the
        # facet-centre visibility test and the missing neighbour patches leave a few percent.
        crater = create_shape_crater(0.4, 0.1; Nx=8, Ny=8)
        shape = load_shape_obj(path_obj; as_hierarchical=true)
        add_roughness_models!(shape, crater)
        T_iso = 250.0
        T_sub = fill(T_iso, length(crater.faces))
        L_lambert = ε * ATPM.σ_SB * T_iso^4 / π
        i = 1
        n̂ = shape.global_shape.face_normals[i]
        @test roughness_radiance(shape, i, T_sub, ε, n̂) ≈ L_lambert rtol=1e-10
        # 30° off the normal, around the facet
        t̂ = normalize(SVector(1.0, 0.0, 0.0) - (SVector(1.0, 0.0, 0.0) ⋅ n̂) * n̂)
        for ϕ in range(0, 2π; length=7)[1:end-1]
            d̂ = cosd(30) * n̂ + sind(30) * (cos(ϕ) * t̂ + sin(ϕ) * (n̂ × t̂))
            @test roughness_radiance(shape, i, T_sub, ε, d̂) ≈ L_lambert rtol=0.05
        end
        @test isnan(roughness_radiance(shape, i, T_sub, ε, -n̂))
        @test_throws ArgumentError roughness_radiance(shape, i, T_sub[1:end-1], ε, n̂)
    end

    @testset "sunlit crater beams towards the Sun" begin
        # The sunlit wall of a crater is hotter and faces the Sun. Seen from the Sun's
        # direction the hot wall is in view; seen from the mirrored direction (same
        # elevation, opposite azimuth) the cold wall is. So T_b(sun) > T_b(anti-sun).
        shape, problem, sol = run(create_shape_crater(0.4, 0.1; Nx=8, Ny=8))
        d̂_sun = SVector(1.0, 0.0, 0.0)
        T_b_sun = brightness_temperature(problem, sol, 1, d̂_sun)
        n_checked = 0
        for i in eachindex(T_b_sun)
            n̂ = shape.global_shape.face_normals[i]
            cosθ = n̂ ⋅ d̂_sun
            0.3 < cosθ < 0.8 || continue           # obliquely lit: beaming is strongest
            d̂_anti = 2 * cosθ * n̂ - d̂_sun          # mirror of d̂_sun about the normal
            T_b_anti = brightness_temperature(problem, sol, 1, d̂_anti)[i]
            T_b_smooth = ε^(1/4) * sol.surface_temperature[i, 1]
            @test T_b_sun[i] > T_b_anti
            @test T_b_sun[i] > T_b_smooth        # beaming raises the apparent temperature
            n_checked += 1
        end
        @test n_checked > 0
    end

    @testset "Planck inversion round-trips" begin
        for T in (150.0, 250.0, 400.0), λ in (5e-6, 10e-6, 20e-6)
            L = ATPM.blackbody_radiance(λ, T) / π
            @test ATPM._brightness_temperature(L, λ) ≈ T rtol=1e-10
        end
        @test ATPM._brightness_temperature(ATPM.σ_SB * 300.0^4 / π, nothing) ≈ 300.0
        # Degenerate radiances: zero is a 0 K blackbody, a negative value is not a radiance
        @test ATPM._brightness_temperature(0.0, 10e-6) == 0.0
        @test isnan(ATPM._brightness_temperature(-1.0, 10e-6))
    end

    @testset "plain ShapeModel and facets without a roughness model are Lambertian" begin
        d̂ = normalize(SVector(1.0, 0.0, 0.0))

        # A plain ShapeModel: every facet is a smooth Lambertian surface
        shape_plain = load_shape_obj(path_obj)
        problem_plain = SingleAsteroidThermoPhysicalProblem(shape_plain, thermo_params, grid_params;
            with_self_shadowing=false, with_self_heating=false)
        output_plain = SingleAsteroidOutputSpec(output_times, Int[]; save_subsurface_temperature=false)
        sol_plain = solve(problem_plain, CrankNicolson(); ephem, output=output_plain,
            initial_temperature=200.0, show_progress=false)
        L_plain = directional_radiance(problem_plain, sol_plain, 1, d̂)
        for i in eachindex(L_plain)
            cosθ = shape_plain.face_normals[i] ⋅ d̂
            if cosθ > 0
                @test L_plain[i] ≈ ε * ATPM.σ_SB * sol_plain.surface_temperature[i, 1]^4 / π
            else
                @test isnan(L_plain[i])
            end
        end

        # A crater on facet 1 only: facet 1 is rough, the others fall back to the smooth
        # Lambertian radiance of their own surface temperature
        shape_part = load_shape_obj(path_obj; as_hierarchical=true)
        add_roughness_models!(shape_part, create_shape_crater(0.4, 0.1; Nx=4, Ny=4), 1)
        problem_part = SingleAsteroidThermoPhysicalProblem(shape_part, thermo_params, grid_params;
            with_self_shadowing=false, with_self_heating=false)
        output_part = SingleAsteroidOutputSpec(output_times, Int[]; save_subsurface_temperature=false,
            roughness_face_ids=[1], save_roughness_surface_temperature=true)
        sol_part = solve(problem_part, CrankNicolson(); ephem, output=output_part,
            initial_temperature=200.0, show_progress=false)
        L_part = directional_radiance(problem_part, sol_part, 1, d̂)
        @test isequal(L_part[2:end], L_plain[2:end])   # same smooth solution, same smooth radiance (NaN-aware)
        @test !isnan(L_part[1]) && L_part[1] != L_plain[1]
    end

    @testset "surface temperature must have been recorded" begin
        shape = load_shape_obj(path_obj; as_hierarchical=true)
        add_roughness_models!(shape, create_shape_crater(0.4, 0.1; Nx=4, Ny=4))
        problem = SingleAsteroidThermoPhysicalProblem(shape, thermo_params, grid_params;
            with_self_shadowing=false, with_self_heating=false)
        output = SingleAsteroidOutputSpec(output_times, Int[];
            save_surface_temperature=false, save_subsurface_temperature=false)
        sol = solve(problem, CrankNicolson(); ephem, output, initial_temperature=200.0, show_progress=false)
        @test_throws ArgumentError directional_radiance(problem, sol, 1, SVector(1.0, 0.0, 0.0))
    end
end
