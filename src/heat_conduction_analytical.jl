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
