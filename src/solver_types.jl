#=
solver_types.jl

Algorithm types for thermophysical model solving, and
cache types for internal heat conduction computations.
=#

# ╔═══════════════════════════════════════════════════════════════════╗
# ║         Algorithm types for thermophysical model solving          ║
# ╚═══════════════════════════════════════════════════════════════════╝

"""
Abstract type for thermophysical model solving algorithms.

Concrete subtypes are passed as the second argument to `solve`:
```julia
solve(problem, CrankNicolson(); kwargs...)
```
"""
abstract type AbstractThermoPhysicalAlgorithm end


"""
    ExplicitEuler()

Explicit (forward) Euler method for heat conduction:
- First order in time
- Conditionally stable (requires Fourier number λ < 0.5)
"""
struct ExplicitEuler <: AbstractThermoPhysicalAlgorithm end


"""
    ImplicitEuler()

Implicit (backward) Euler method for heat conduction:
- First order in time
- Unconditionally stable
"""
struct ImplicitEuler <: AbstractThermoPhysicalAlgorithm end


"""
    CrankNicolson()

Crank-Nicolson method for heat conduction:
- Second order in time
- Unconditionally stable

# References
- https://en.wikipedia.org/wiki/Crank–Nicolson_method
"""
struct CrankNicolson <: AbstractThermoPhysicalAlgorithm end


# ╔═══════════════════════════════════════════════════════════════════╗
# ║       Cache types for internal heat conduction computations       ║
# ╚═══════════════════════════════════════════════════════════════════╝

"""
Abstract type for internal pre-allocated caches used in heat conduction computations.
These types are not part of the public API.
"""
abstract type HeatConductionCache end


"""
Internal cache for the explicit (forward) Euler method.
Holds a pre-allocated vector for the temperature at the next time step.
"""
struct ExplicitEulerCache <: HeatConductionCache
    x::Vector{Float64}
end

ExplicitEulerCache(grid_params::GridParams) = ExplicitEulerCache(grid_params.n_depth)
ExplicitEulerCache(N::Integer) = ExplicitEulerCache(zeros(N))


"""
Internal cache for the implicit (backward) Euler method.
Holds pre-allocated vectors for the tridiagonal matrix algorithm.
"""
struct ImplicitEulerCache <: HeatConductionCache
    a::Vector{Float64}
    b::Vector{Float64}
    c::Vector{Float64}
    d::Vector{Float64}
    x::Vector{Float64}
end

ImplicitEulerCache(grid_params::GridParams) = ImplicitEulerCache(grid_params.n_depth)
ImplicitEulerCache(N::Integer) = ImplicitEulerCache(zeros(N), zeros(N), zeros(N), zeros(N), zeros(N))


"""
Internal cache for the Crank-Nicolson method.
Holds pre-allocated vectors for the tridiagonal matrix algorithm.

# References
- https://en.wikipedia.org/wiki/Crank–Nicolson_method
"""
struct CrankNicolsonCache <: HeatConductionCache
    a::Vector{Float64}
    b::Vector{Float64}
    c::Vector{Float64}
    d::Vector{Float64}
    x::Vector{Float64}
end

CrankNicolsonCache(grid_params::GridParams) = CrankNicolsonCache(grid_params.n_depth)
CrankNicolsonCache(N::Integer) = CrankNicolsonCache(zeros(N), zeros(N), zeros(N), zeros(N), zeros(N))


# ╔═══════════════════════════════════════════════════════════════════╗
# ║    Types of boundary conditions for heat conduction equations    ║
# ╚═══════════════════════════════════════════════════════════════════╝

"""
Abstract type of a boundary condition for a heat conduction equation
"""
abstract type AbstractBoundaryCondition end


"""
Singleton type of radiation boundary condition
"""
struct RadiationBoundaryCondition <: AbstractBoundaryCondition end


"""
Singleton type of insulation boundary condition
"""
struct InsulationBoundaryCondition <: AbstractBoundaryCondition end


"""
Type of isothermal boundary condition
"""
struct IsothermalBoundaryCondition <: AbstractBoundaryCondition
    T_iso::Float64
end
