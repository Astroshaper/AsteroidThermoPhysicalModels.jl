#=
heat_conduction_analytical.jl

Analytical solutions for the 1D heat conduction equation.
This file provides exact solutions for specific boundary and initial conditions,
useful for validating numerical solvers and understanding thermal behavior.
=#

"""
    analytical_solution_isothermal(x, t, L, α; n_max=100) -> T

Calculate the analytical solution of the 1D heat equation with isothermal boundary conditions.
- Equation: ∂T/∂t = α ∂²T/∂x²
- Domain: 0 ≤ x ≤ L
- Boundary conditions: T(0,t) = T(L,t) = 0
- Initial condition: T₀(x) = x < 0.5L ? 2x/L : 2(1 - x/L)  # Triangular profile as follows:
        T₀
        ^ 
      1 |   ・
        |  ・ ・
        | ・   ・
        |・     ・
      0 +---+---+--> x
        0  L/2  L

The solution is given by the Fourier series:
    T(x, t) = Σ Bₙ * sin(nπx/L) * exp(-αn²π²t/L²)
where Bₙ = (2/L) * ∫₀^L T₀(ξ) * sin(nπξ/L) dξ
For the triangular initial condition, the coefficients can be calculated analytically:
    Bₙ = (8/n²π²) * sin(nπ/2)
only for odd `n`. The sum of even-`n` terms is zero due to symmetry.

# Arguments
- `x`     : Position [m]
- `t`     : Time [s]
- `L`     : Length of the domain [m]
- `α`     : Thermal diffusivity [m²/s]
- `n_max` : Number of terms in the Fourier series

# Returns
- `T`     : Temperature [K]
"""
function analytical_solution_isothermal(x, t, L, α; n_max=100)
    T = 0.0
    for n in 1:2:n_max  # Only odd terms contribute for this initial condition
        Bₙ = 8/ (n * π)^2 * sin(n*π/2)
        T += Bₙ * sin(n*π*x/L) * exp(-α*(n*π/L)^2*t)
    end
    return T
end


"""
    analytical_solution_periodic(z, L, P, k, ρ, Cₚ, ε, T₀, F₁) -> T̃

Complex amplitude of the periodic (linearised) solution of the 1D heat equation in a slab of
thickness `L` with a radiation boundary condition at the surface and an insulated bottom.
- Equation: ∂T/∂t = α ∂²T/∂z², α = k/(ρCₚ), 0 ≤ z ≤ L
- Surface (z = 0): F_abs(t) = εσT⁴ - k ∂T/∂z, with F_abs(t) = εσT₀⁴ + Re[F₁ exp(iωt)]
- Bottom (z = L): ∂T/∂z = 0
- ω = 2π/P

The radiative term is linearised around `T₀`, εσT⁴ ≈ εσT₀⁴ + 4εσT₀³ (T - T₀), which is
accurate to O((ΔT/T₀)²). The temperature is then

    T(z, t) = T₀ + Re[T̃(z) exp(iωt)],    T̃(z) = A cosh(q(L - z)) / cosh(qL)

with q = (1 + i)/l_s, l_s = √(2α/ω) the e-folding skin depth, and

    A = ΔT_eq / (1 + (1 + i) β tanh(qL)),    ΔT_eq = F₁ / (4εσT₀³),    β = Γ√ω / (4√2 εσT₀³)

where Γ = √(kρCₚ) is the thermal inertia. `|T̃(z)|` is the amplitude of the temperature
oscillation at depth `z` and `-arg(T̃(z))` its phase lag behind the forcing. For `L → ∞`
the amplitude decays as `exp(-z/l_s)` and the lag grows as `z/l_s`; the surface lag never
exceeds 45°.

# Arguments
- `z`  : Depth [m]
- `L`  : Thickness of the slab [m]
- `P`  : Period of the forcing [s]
- `k`  : Thermal conductivity [W/m/K]
- `ρ`  : Density [kg/m³]
- `Cₚ` : Heat capacity [J/kg/K]
- `ε`  : Emissivity [-]
- `T₀` : Mean temperature, εσT₀⁴ = ⟨F_abs⟩ [K]
- `F₁` : Complex amplitude of the absorbed flux [W/m²] (real for a forcing in phase with cos(ωt))

# Returns
- `T̃` : Complex amplitude of the temperature at depth `z` [K]

# References
- Carslaw & Jaeger (1959), *Conduction of Heat in Solids*, 2nd ed., §3.6
- Spencer et al. (1989), Icarus 78, 337 (the thermal parameter Θ = 4√2 β)
"""
function analytical_solution_periodic(z, L, P, k, ρ, Cₚ, ε, T₀, F₁)
    ω  = 2π / P
    α  = thermal_diffusivity(k, ρ, Cₚ)
    Γ  = thermal_inertia(k, ρ, Cₚ)
    εσ = ε * σ_SB

    l_s = √(2α / ω)          # e-folding skin depth
    q   = (1 + im) / l_s     # complex wavenumber

    ΔT_eq = F₁ / (4εσ * T₀^3)
    β     = Γ * √ω / (4√2 * εσ * T₀^3)
    A     = ΔT_eq / (1 + (1 + im) * β * tanh(q * L))

    return A * cosh(q * (L - z)) / cosh(q * L)
end
