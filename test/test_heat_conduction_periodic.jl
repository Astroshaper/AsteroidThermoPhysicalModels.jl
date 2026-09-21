#=
test_heat_conduction_periodic.jl

Verification of the 1D heat conduction solvers under the radiation boundary condition
against the periodic analytic solution (linearised radiation BC, insulated bottom):
- Amplitude and phase lag of the surface temperature
- Decay and phase lag with depth
- Error decreases when the grid is refined
- Time-averaged energy balance εσ⟨T_s⁴⟩ = ⟨F_abs⟩ in the periodic state

The thresholds are set to what the current solvers achieve and document their accuracy;
tighten them when the surface boundary condition is improved (#247).
=#

@testset "heat_conduction_periodic" begin
    msg = """\n
    ┏━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━┓
    |             Test: heat_conduction_periodic             |
    ┗━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━┛
    """
    println(msg)

    ATPM = AsteroidThermoPhysicalModels

    ## --- Shape model (a single facet) ---
    shape = load_shape_obj(joinpath("shape", "single_face.obj"))

    ## --- Thermal properties ---
    k  = 0.1     # Thermal conductivity [W/m/K]
    ρ  = 1270.0  # Density [kg/m³]
    Cₚ = 600.0   # Heat capacity [J/kg/K]
    ε  = 1.0     # Emissivity [-]
    thermo_params = ThermoParams(k, ρ, Cₚ, 0.0, 0.0, ε)

    ## --- Forcing: F_abs(t) = F̄ (1 + a cos ωt) ---
    P  = 7.633 * 3600           # Period [s]
    ω  = 2π / P
    T₀ = 300.0                  # Mean temperature [K]
    F̄  = ε * ATPM.σ_SB * T₀^4   # Mean absorbed flux [W/m²]
    a  = 0.08                   # Relative forcing amplitude; ΔT/T₀ ≈ 0.013 keeps the linearisation error ≲ 1e-4
    F₁ = F̄ * a

    ## --- Domain: ζ = L/l_s skin depths deep ---
    α   = ATPM.thermal_diffusivity(k, ρ, Cₚ)
    l_s = √(2α / ω)             # e-folding skin depth [m]
    ζ   = 5.0
    L   = ζ * l_s

    n_spinup = 20  # Periods before sampling; the slowest mode relaxes as exp(-t / 0.129 ζ² P)

    wrap_deg(x) = mod(x + 180, 360) - 180  # Wrap an angle difference into (-180, 180]

    """
        run_periodic(algorithm, step!, N_z, N_t) -> (zs, c₁, εσT⁴_mean, F_mean)

    Spin up a single column for `n_spinup` periods and, over one more period, accumulate
    the first Fourier coefficient `c₁` of the temperature at every depth node, the mean
    of εσT_s⁴ and the mean absorbed flux.
    - `N_z` : Nodes per skin depth
    - `N_t` : Time steps per period
    """
    function run_periodic(algorithm, step!, N_z, N_t)
        n_depth = round(Int, ζ * N_z) + 1
        Δt = P / N_t
        grid_params = GridParams(; z_max = L, n_depth)
        problem = SingleAsteroidThermoPhysicalProblem(shape, thermo_params, grid_params;
            with_self_shadowing      = false,
            with_self_heating        = false,
            upper_boundary_condition = RadiationBoundaryCondition(),
            lower_boundary_condition = InsulationBoundaryCondition(),
        )
        state = ATPM._build_single_state(problem, algorithm)
        state.temperature .= T₀

        for n in 0:(n_spinup * N_t - 1)
            state.flux_sun[1] = F̄ * (1 + a * cos(ω * n * Δt))
            step!(state, Δt)
        end

        c₁ = zeros(ComplexF64, n_depth)
        εσT⁴_mean = 0.0
        F_mean = 0.0
        for n in 0:(N_t - 1)
            t = (n_spinup * N_t + n) * Δt
            state.flux_sun[1] = F̄ * (1 + a * cos(ω * t))
            step!(state, Δt)
            c₁ .+= (2 / N_t) .* state.temperature[:, 1] .* exp(-im * ω * (t + Δt))  # T is now at t + Δt
            εσT⁴_mean += ε * ATPM.σ_SB * state.temperature[1, 1]^4 / N_t
            F_mean    += state.flux_sun[1] / N_t
        end

        zs = grid_params.Δz .* (0:n_depth-1)
        return zs, c₁, εσT⁴_mean, F_mean
    end

    solvers = (
        ("Explicit Euler", ExplicitEuler(), ATPM.explicit_euler!),
        ("Implicit Euler", ImplicitEuler(), ATPM.implicit_euler!),
        ("Crank-Nicolson", CrankNicolson(), ATPM.crank_nicolson!),
    )

    T̃(z) = ATPM.analytical_solution_periodic(z, L, P, k, ρ, Cₚ, ε, T₀, F₁)
    A = T̃(0.0)
    println("Analytic surface temperature: amplitude = $(abs(A)) K, phase lag = $(-rad2deg(angle(A)))°")

    ## Surface amplitude errors at N_z = 10 and 20, kept for the convergence check
    amp_err = Dict{String, Vector{Float64}}()

    for (name, algorithm, step!) in solvers
        @testset "$name vs. periodic analytic solution" begin
            N_z, N_t = 20, 4000  # λ = αΔt/Δz² = π N_z² / N_t = 0.31, stable for the explicit solver
            zs, c₁, εσT⁴_mean, F_mean = run_periodic(algorithm, step!, N_z, N_t)

            ## Surface amplitude and phase lag
            δ_amp = abs(c₁[1]) / abs(A) - 1
            δ_lag = wrap_deg(rad2deg(angle(A) - angle(c₁[1])))
            println("$name: surface amplitude error = $(δ_amp), phase lag error = $(δ_lag)°")
            @test abs(δ_amp) < 1e-2
            @test abs(δ_lag) < 1.0

            ## Decay and phase lag with depth, relative to the surface
            for x in (1.0, 2.0, 3.0)
                i = round(Int, x * N_z) + 1
                ratio_num = abs(c₁[i]) / abs(c₁[1])
                ratio_ana = abs(T̃(zs[i])) / abs(A)
                lag_num = wrap_deg(rad2deg(angle(c₁[1]) - angle(c₁[i])))
                lag_ana = wrap_deg(rad2deg(angle(A) - angle(T̃(zs[i]))))
                @test abs(ratio_num / ratio_ana - 1) < 5e-3
                @test abs(wrap_deg(lag_num - lag_ana)) < 0.2
            end

            ## Time-averaged energy balance in the periodic state
            @test abs(εσT⁴_mean / F_mean - 1) < 1e-5

            amp_err[name] = [abs(δ_amp)]
        end
    end

    @testset "Error decreases with grid refinement" begin
        for (name, algorithm, step!) in solvers
            _, c₁, _, _ = run_periodic(algorithm, step!, 10, 1000)  # same λ as above
            pushfirst!(amp_err[name], abs(abs(c₁[1]) / abs(A) - 1))
            ratio = amp_err[name][1] / amp_err[name][2]
            println("$name: surface amplitude error N_z = 10 → 20: $(amp_err[name][1]) → $(amp_err[name][2]) (ratio $ratio)")
            @test ratio > 1.5
        end
    end
end
