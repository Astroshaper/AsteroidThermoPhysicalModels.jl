#=
tpm_problem.jl

Problem types for thermophysical model solving.
A problem type encapsulates everything needed to define "what to solve":
shape model, thermophysical parameters, boundary conditions, and modeling flags.
=#

"""
Abstract type for thermophysical problem definitions.
"""
abstract type AbstractThermoPhysicalProblem end


# ╔═══════════════════════════════════════════════════════════════════╗
# ║              Single-asteroid thermophysical problem               ║
# ╚═══════════════════════════════════════════════════════════════════╝

"""
    struct SingleAsteroidThermoPhysicalProblem

Defines the thermophysical problem for a single asteroid.
Encapsulates all information needed to describe the physical problem,
separate from the numerical method used to solve it.

# Fields
- `shape`                    : Shape model of the asteroid (`ShapeModel` or `HierarchicalShapeModel`)
- `thermo_params`            : Thermophysical material parameters (all vectors expanded to length `n_face`)
- `grid_params`              : Numerical grid settings
- `with_self_shadowing`      : Whether to include self-shadowing
- `with_self_heating`        : Whether to include self-heating (re-absorption of thermal emission from other faces)
- `upper_boundary_condition` : Boundary condition at the surface (upper boundary)
- `lower_boundary_condition` : Boundary condition at depth (lower boundary)

# Usage
```julia
thermo_params = ThermoParams(k, ρ, Cₚ, R_vis, R_ir, ε)
grid_params   = GridParams(z_max, n_depth, Δz)
problem = SingleAsteroidThermoPhysicalProblem(shape, thermo_params, grid_params;
    with_self_shadowing = true,
    with_self_heating   = true,
    upper_boundary_condition = RadiationBoundaryCondition(),
    lower_boundary_condition = InsulationBoundaryCondition(),
)
solution = solve(problem, CrankNicolson(); ephem = ephem, initial_temperature = 200.0)
```
"""
struct SingleAsteroidThermoPhysicalProblem{
    Sh <: AbstractShapeModel,
    BU <: AbstractBoundaryCondition,
    BL <: AbstractBoundaryCondition,
} <: AbstractThermoPhysicalProblem
    shape                    ::Sh
    thermo_params            ::ThermoParams
    grid_params              ::GridParams
    with_self_shadowing      ::Bool
    with_self_heating        ::Bool
    upper_boundary_condition ::BU
    lower_boundary_condition ::BL
end


"""
    _expand_thermo_params(thermo_params::ThermoParams, n_face::Int) -> ThermoParams

Expand a `ThermoParams` with length-1 vectors to length `n_face` (uniform surface),
or validate that all vectors already have length `n_face` (non-uniform surface).
"""
function _expand_thermo_params(thermo_params::ThermoParams, n_face::Int)
    n = length(thermo_params.conductivity)
    if n == 1
        return ThermoParams(
            fill(thermo_params.conductivity[1],    n_face),
            fill(thermo_params.density[1],         n_face),
            fill(thermo_params.heat_capacity[1],   n_face),
            fill(thermo_params.reflectance_vis[1], n_face),
            fill(thermo_params.reflectance_ir[1],  n_face),
            fill(thermo_params.emissivity[1],      n_face),
        )
    elseif n == n_face
        return thermo_params
    else
        throw(ArgumentError("ThermoParams vector length ($n) must be 1 (uniform) or n_face ($n_face)."))
    end
end


"""
    SingleAsteroidThermoPhysicalProblem(shape, thermo_params, grid_params; kwargs...) -> problem

Construct a thermophysical problem for a single asteroid.

# Arguments
- `shape`         : Shape model of the asteroid
- `thermo_params` : Thermophysical material parameters (`ThermoParams`)
- `grid_params`   : Numerical grid settings (`GridParams`)

# Keyword Arguments
- `with_self_shadowing = true` : Whether to include self-shadowing
- `with_self_heating   = true` : Whether to include self-heating
- `upper_boundary_condition = RadiationBoundaryCondition()`  : Boundary condition at the surface
- `lower_boundary_condition = InsulationBoundaryCondition()` : Boundary condition at depth

# Notes
- If `with_self_shadowing = true` and `face_visibility_graph` is not yet built, it is built automatically.
  To avoid this, pre-build with `build_face_visibility_graph!` or pass `with_face_visibility=true`
  when loading the shape.
- `ThermoParams` with length-1 vectors is expanded to `n_face` at construction time.
"""
function SingleAsteroidThermoPhysicalProblem(shape, thermo_params::ThermoParams, grid_params::GridParams;
    with_self_shadowing ::Bool = true,
    with_self_heating   ::Bool = true,
    upper_boundary_condition   = RadiationBoundaryCondition(),
    lower_boundary_condition   = InsulationBoundaryCondition(),
)
    if with_self_shadowing && isnothing(shape.face_visibility_graph)
        @info "Building face_visibility_graph for self-shadowing..."
        build_face_visibility_graph!(shape)
    end

    n_face = length(shape.faces)
    thermo_params_expanded = _expand_thermo_params(thermo_params, n_face)

    SingleAsteroidThermoPhysicalProblem(shape, thermo_params_expanded, grid_params, with_self_shadowing, with_self_heating, upper_boundary_condition, lower_boundary_condition)
end


# ╔═══════════════════════════════════════════════════════════════════╗
# ║               Binary-asteroid thermophysical problem              ║
# ╚═══════════════════════════════════════════════════════════════════╝

"""
    struct BinaryAsteroidThermoPhysicalProblem

Defines the thermophysical problem for a binary asteroid system.

# Fields
- `primary`               : Problem definition for the primary body
- `secondary`             : Problem definition for the secondary body
- `with_mutual_shadowing` : Whether to include mutual shadowing (eclipses)
- `with_mutual_heating`   : Whether to include mutual heating

# Usage
```julia
prob1 = SingleAsteroidThermoPhysicalProblem(shape1, thermo_params1; ...)
prob2 = SingleAsteroidThermoPhysicalProblem(shape2, thermo_params2; ...)
problem = BinaryAsteroidThermoPhysicalProblem(prob1, prob2;
    with_mutual_shadowing = true,
    with_mutual_heating   = true,
)
solution = solve(problem, CrankNicolson(); ephem = ephem, T₀ = 200.0)
```
"""
struct BinaryAsteroidThermoPhysicalProblem{
    P1 <: SingleAsteroidThermoPhysicalProblem,
    P2 <: SingleAsteroidThermoPhysicalProblem,
} <: AbstractThermoPhysicalProblem
    primary               ::P1
    secondary             ::P2
    with_mutual_shadowing ::Bool
    with_mutual_heating   ::Bool
end


"""
    BinaryAsteroidThermoPhysicalProblem(primary, secondary; kwargs...) -> problem

Construct a thermophysical problem for a binary asteroid system.

# Arguments
- `primary`   : Problem definition for the primary body
- `secondary` : Problem definition for the secondary body

# Keyword Arguments
- `with_mutual_shadowing = true` : Whether to include mutual shadowing (eclipses)
- `with_mutual_heating   = true` : Whether to include mutual heating

# Notes
- If `with_mutual_shadowing = true` and BVH is not yet built for either shape, it is built automatically.
  To avoid this, pre-build with `build_bvh!` or pass `with_bvh=true` when loading shapes.
"""
function BinaryAsteroidThermoPhysicalProblem(
    primary   ::SingleAsteroidThermoPhysicalProblem,
    secondary ::SingleAsteroidThermoPhysicalProblem;
    with_mutual_shadowing ::Bool = true,
    with_mutual_heating   ::Bool = true,
)
    if with_mutual_shadowing
        if isnothing(primary.shape.bvh)
            @info "Building BVH for primary shape..."
            build_bvh!(primary.shape)
        end
        if isnothing(secondary.shape.bvh)
            @info "Building BVH for secondary shape..."
            build_bvh!(secondary.shape)
        end
    end

    BinaryAsteroidThermoPhysicalProblem(primary, secondary, with_mutual_shadowing, with_mutual_heating)
end


_to_pair(x)        = (x, x)   # single value → apply to both bodies
_to_pair(x::Tuple) = x        # Tuple → use as-is (element per body)


"""
    BinaryAsteroidThermoPhysicalProblem(shape, thermo_params, grid_params; kwargs...) -> problem

Convenience constructor: build both single-body problems from raw shapes and parameters.

# Arguments
- `shape`         : Tuple `(shape1, shape2)` of shape models
- `thermo_params` : `ThermoParams` applied to both bodies, or Tuple `(thermo_params1, thermo_params2)` for per-body settings
- `grid_params`   : `GridParams` applied to both bodies, or Tuple `(grid_params1, grid_params2)` for per-body settings

# Keyword Arguments
Single-body kwargs (applied identically to both bodies):
- `with_self_shadowing      = true`
- `with_self_heating        = true`
- `upper_boundary_condition = RadiationBoundaryCondition()`
- `lower_boundary_condition = InsulationBoundaryCondition()`

Binary-system kwargs:
- `with_mutual_shadowing = true`
- `with_mutual_heating   = true`

# Notes
For per-body control of single-body kwargs, use `SingleAsteroidThermoPhysicalProblem`
separately and pass the results to `BinaryAsteroidThermoPhysicalProblem(primary, secondary; ...)`.

# Example
```julia
grid_params = GridParams(z_max, Δz, n_depth)
problem = BinaryAsteroidThermoPhysicalProblem(
    (shape1, shape2),
    (thermo_params1, thermo_params2),
    grid_params;
    with_mutual_shadowing = true,
    with_mutual_heating   = true,
)
```
"""
function BinaryAsteroidThermoPhysicalProblem(
    shape         ::Tuple,
    thermo_params ::Union{Tuple, ThermoParams},
    grid_params   ::Union{Tuple, GridParams};
    with_self_shadowing::Bool   = true,
    with_self_heating::Bool     = true,
    upper_boundary_condition    = RadiationBoundaryCondition(),
    lower_boundary_condition    = InsulationBoundaryCondition(),
    with_mutual_shadowing::Bool = true,
    with_mutual_heating::Bool   = true,
)
    thermo_params1, thermo_params2 = _to_pair(thermo_params)
    grid_params1, grid_params2 = _to_pair(grid_params)

    prob1 = SingleAsteroidThermoPhysicalProblem(shape[1], thermo_params1, grid_params1;
        with_self_shadowing,
        with_self_heating,
        upper_boundary_condition,
        lower_boundary_condition,
    )

    prob2 = SingleAsteroidThermoPhysicalProblem(shape[2], thermo_params2, grid_params2;
        with_self_shadowing,
        with_self_heating,
        upper_boundary_condition,
        lower_boundary_condition,
    )

    BinaryAsteroidThermoPhysicalProblem(prob1, prob2;
        with_mutual_shadowing,
        with_mutual_heating,
    )
end
