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


"""
    struct SingleAsteroidThermoPhysicalProblem

Defines the thermophysical problem for a single asteroid.
Encapsulates all information needed to describe the physical problem,
separate from the numerical method used to solve it.

# Fields
- `shape`                    : Shape model of the asteroid (`ShapeModel` or `HierarchicalShapeModel`)
- `thermo_params`            : Thermophysical parameters
- `with_self_shadowing`      : Whether to include self-shadowing
- `with_self_heating`        : Whether to include self-heating (re-absorption of thermal emission from other faces)
- `upper_boundary_condition` : Boundary condition at the surface (upper boundary)
- `lower_boundary_condition` : Boundary condition at depth (lower boundary)

# Usage
```julia
problem = SingleAsteroidThermoPhysicalProblem(shape, thermo_params;
    with_self_shadowing = true,
    with_self_heating   = true,
    upper_boundary_condition = RadiationBoundaryCondition(),
    lower_boundary_condition = InsulationBoundaryCondition(),
)
solution = solve(problem, CrankNicolson(); ephem = ephem, T₀ = 200.0)
```
"""
struct SingleAsteroidThermoPhysicalProblem{
    Sh <: AbstractShapeModel,
    P  <: AbstractThermoParams,
    BU <: AbstractBoundaryCondition,
    BL <: AbstractBoundaryCondition,
} <: AbstractThermoPhysicalProblem
    shape                    ::Sh
    thermo_params            ::P
    with_self_shadowing      ::Bool
    with_self_heating        ::Bool
    upper_boundary_condition ::BU
    lower_boundary_condition ::BL
end


"""
    SingleAsteroidThermoPhysicalProblem(shape, thermo_params; kwargs...) -> problem

Construct a thermophysical problem for a single asteroid.

# Arguments
- `shape`         : Shape model of the asteroid
- `thermo_params` : Thermophysical parameters

# Keyword Arguments
- `with_self_shadowing = true` : Whether to include self-shadowing
- `with_self_heating   = true` : Whether to include self-heating
- `upper_boundary_condition = RadiationBoundaryCondition()`  : Boundary condition at the surface
- `lower_boundary_condition = InsulationBoundaryCondition()` : Boundary condition at depth

# Notes
- If `with_self_shadowing = true` and `face_visibility_graph` is not yet built, it is built automatically.
  To avoid this, pre-build with `build_face_visibility_graph!` or pass `with_face_visibility=true`
  when loading the shape.
- Thermophysical parameters are broadcast to all faces via `broadcast_thermo_params!`.
"""
function SingleAsteroidThermoPhysicalProblem(shape, thermo_params;
    with_self_shadowing = true,
    with_self_heating   = true,
    upper_boundary_condition = RadiationBoundaryCondition(),
    lower_boundary_condition = InsulationBoundaryCondition(),
)
    if with_self_shadowing && isnothing(shape.face_visibility_graph)
        @info "Building face_visibility_graph for self-shadowing..."
        build_face_visibility_graph!(shape)
    end

    broadcast_thermo_params!(thermo_params, shape)

    SingleAsteroidThermoPhysicalProblem(shape, thermo_params, with_self_shadowing, with_self_heating, upper_boundary_condition, lower_boundary_condition)
end


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
prob_primary   = SingleAsteroidThermoPhysicalProblem(shape1, params1; ...)
prob_secondary = SingleAsteroidThermoPhysicalProblem(shape2, params2; ...)
problem = BinaryAsteroidThermoPhysicalProblem(prob_primary, prob_secondary;
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
function BinaryAsteroidThermoPhysicalProblem(primary, secondary;
    with_mutual_shadowing = true,
    with_mutual_heating   = true,
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
