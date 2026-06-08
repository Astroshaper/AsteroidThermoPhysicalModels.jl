#=
solver_types.jl

Types for heat conduction solvers and boundary conditions.
=#

# ╔═══════════════════════════════════════════════════════════════════╗
# ║          Types of solvers for heat conduction equations           ║
# ╚═══════════════════════════════════════════════════════════════════╝

"""
Abstract type of a solver for a heat conduction equation
"""
abstract type HeatConductionSolver end


"""
Type of the explicit (forward) Euler method:
- Explicit in time
- First order in time

The `ExplicitEulerSolver` type includes a vector for the temperature at the next time step.
"""
struct ExplicitEulerSolver <: HeatConductionSolver
    x::Vector{Float64}  # Temperature vector for the next time step
end

ExplicitEulerSolver(thermo_params::AbstractThermoParams) = ExplicitEulerSolver(thermo_params.n_depth)
ExplicitEulerSolver(N::Integer) = ExplicitEulerSolver(zeros(N))


"""
Type of the implicit (backward) Euler method:
- Implicit in time (Unconditionally stable in the heat conduction equation)
- First order in time

The `ImplicitEulerSolver` type has vectors for the tridiagonal matrix algorithm.
"""
struct ImplicitEulerSolver <: HeatConductionSolver
    a::Vector{Float64}
    b::Vector{Float64}
    c::Vector{Float64}
    d::Vector{Float64}
    x::Vector{Float64}
end

ImplicitEulerSolver(thermo_params::AbstractThermoParams) = ImplicitEulerSolver(thermo_params.n_depth)
ImplicitEulerSolver(N::Integer) = ImplicitEulerSolver(zeros(N), zeros(N), zeros(N), zeros(N), zeros(N))


"""
Type of the Crank-Nicolson method:
- Implicit in time (Unconditionally stable in the heat conduction equation)
- Second order in time

The `CrankNicolsonSolver` type has vectors for the tridiagonal matrix algorithm.

# References
- https://en.wikipedia.org/wiki/Crank–Nicolson_method
"""
struct CrankNicolsonSolver <: HeatConductionSolver
    a::Vector{Float64}
    b::Vector{Float64}
    c::Vector{Float64}
    d::Vector{Float64}
    x::Vector{Float64}
end

CrankNicolsonSolver(thermo_params::AbstractThermoParams) = CrankNicolsonSolver(thermo_params.n_depth)
CrankNicolsonSolver(N::Integer) = CrankNicolsonSolver(zeros(N), zeros(N), zeros(N), zeros(N), zeros(N))


# ╔═══════════════════════════════════════════════════════════════════╗
# ║    Types of boundary conditions for heat conduction equations    ║
# ╚═══════════════════════════════════════════════════════════════════╝

"""
Abstract type of a boundary condition for a heat conduction equation
"""
abstract type BoundaryCondition end


"""
Singleton type of radiation boundary condition
"""
struct RadiationBoundaryCondition <: BoundaryCondition end


"""
Singleton type of insulation boundary condition
"""
struct InsulationBoundaryCondition <: BoundaryCondition end


"""
Type of isothermal boundary condition
"""
struct IsothermalBoundaryCondition <: BoundaryCondition
    T_iso::Float64
end
