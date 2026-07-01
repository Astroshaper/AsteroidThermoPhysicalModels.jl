#=
test_thermo_grid_params.jl

Unit tests for ThermoParams and GridParams types:
- Construction from scalars, vectors, and mixed arguments
- Keyword argument constructors
- Automatic broadcast of scalar arguments to vector length
- Validation of vector length consistency
=#

@testset "ThermoParams" begin
    msg = """
    ┏━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━┓
    |               Test: ThermoParams                      |
    ┗━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━┛
    """
    println(msg)

    k     = 0.1
    ρ     = 1500.0
    Cₚ    = 800.0
    R_vis = 0.05
    R_ir  = 0.0
    ε     = 0.9

    @testset "Scalar constructor — uniform surface" begin
        tp = ThermoParams(k, ρ, Cₚ, R_vis, R_ir, ε)

        @test tp isa ThermoParams
        @test length(tp.conductivity)    == 1
        @test tp.conductivity[1]    ≈ k
        @test tp.density[1]         ≈ ρ
        @test tp.heat_capacity[1]   ≈ Cₚ
        @test tp.reflectance_vis[1] ≈ R_vis
        @test tp.reflectance_ir[1]  ≈ R_ir
        @test tp.emissivity[1]      ≈ ε
    end

    @testset "Keyword constructor — uniform surface" begin
        tp = ThermoParams(
            conductivity    = k,
            density         = ρ,
            heat_capacity   = Cₚ,
            reflectance_vis = R_vis,
            reflectance_ir  = R_ir,
            emissivity      = ε,
        )

        @test tp isa ThermoParams
        @test tp.conductivity[1]    ≈ k
        @test tp.reflectance_vis[1] ≈ R_vis
        @test tp.emissivity[1]      ≈ ε
    end

    @testset "Vector constructor — non-uniform surface" begin
        n = 5
        k_vec     = fill(k, n)
        R_vis_vec = [0.04, 0.05, 0.06, 0.07, 0.08]
        tp = ThermoParams(k_vec, fill(ρ, n), fill(Cₚ, n), R_vis_vec, fill(R_ir, n), fill(ε, n))

        @test length(tp.conductivity)    == n
        @test length(tp.reflectance_vis) == n
        @test tp.reflectance_vis ≈ R_vis_vec
    end

    @testset "Mixed scalar/vector constructor" begin
        n = 4
        k_vec = [0.1, 0.2, 0.3, 0.4]
        tp = ThermoParams(k_vec, ρ, Cₚ, R_vis, R_ir, ε)

        @test length(tp.conductivity)  == n
        @test length(tp.density)       == n
        @test length(tp.emissivity)    == n
        @test tp.conductivity         ≈ k_vec
        @test all(tp.density         .≈ ρ)
        @test all(tp.heat_capacity   .≈ Cₚ)
        @test all(tp.emissivity      .≈ ε)
    end

    @testset "Keyword constructor — mixed scalar/vector" begin
        n = 3
        k_vec = [0.1, 0.2, 0.3]
        tp = ThermoParams(
            conductivity    = k_vec,
            density         = ρ,
            heat_capacity   = Cₚ,
            reflectance_vis = R_vis,
            reflectance_ir  = R_ir,
            emissivity      = ε,
        )

        @test length(tp.conductivity) == n
        @test all(tp.density .≈ ρ)
    end

    @testset "Validation — mismatched vector lengths" begin
        k_vec = [0.1, 0.2, 0.3]
        ρ_vec = [1500.0, 1600.0]  # different length
        @test_throws ArgumentError ThermoParams(k_vec, ρ_vec, Cₚ, R_vis, R_ir, ε)
    end
end


@testset "GridParams" begin
    msg = """
    ┏━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━┓
    |               Test: GridParams                        |
    ┗━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━┛
    """
    println(msg)

    z_max   = 0.6
    n_depth = 61

    @testset "Keyword constructor — auto Δz" begin
        gp = GridParams(; z_max, n_depth)

        @test gp isa GridParams
        @test gp.z_max   ≈ z_max
        @test gp.n_depth == n_depth
        @test gp.Δz      ≈ z_max / (n_depth - 1)
    end

    @testset "Positional constructor — explicit Δz" begin
        Δz = z_max / (n_depth - 1)
        gp = GridParams(z_max, n_depth, Δz)

        @test gp.z_max   ≈ z_max
        @test gp.n_depth == n_depth
        @test gp.Δz      ≈ Δz
    end

    @testset "Both constructors are consistent" begin
        gp_kw  = GridParams(; z_max, n_depth)
        gp_pos = GridParams(z_max, n_depth, z_max / (n_depth - 1))

        @test gp_kw.z_max   ≈ gp_pos.z_max
        @test gp_kw.n_depth == gp_pos.n_depth
        @test gp_kw.Δz      ≈ gp_pos.Δz
    end
end
