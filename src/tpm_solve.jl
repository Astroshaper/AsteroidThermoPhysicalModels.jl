#=
tpm_solve.jl

Implementation of the solve interface for thermophysical problems.
Extends CommonSolve.solve for AbstractThermoPhysicalProblem types.

Public API:
    solve(problem, algorithm; ephem, ...) -> solution

Internal dispatch:
    _solve(problem, algorithm, ephem::ConcreteEphemType; ...) -> solution
=#

# Build the appropriate HeatConductionCache from an algorithm type
_build_cache(::ExplicitEuler, grid_params) = ExplicitEulerCache(grid_params)
_build_cache(::ImplicitEuler, grid_params) = ImplicitEulerCache(grid_params)
_build_cache(::CrankNicolson, grid_params) = CrankNicolsonCache(grid_params)


# Internal helper: build SingleAsteroidThermoPhysicalState from a problem + algorithm.
function _build_single_state(
    problem   ::SingleAsteroidThermoPhysicalProblem,
    algorithm ::AbstractThermoPhysicalAlgorithm,
)
    cache   = _build_cache(algorithm, problem.grid_params)
    n_depth = problem.grid_params.n_depth
    n_face  = length(problem.shape.faces)

    SingleAsteroidThermoPhysicalState(
        problem,
        cache,
        zeros(Bool, n_face),
        zeros(n_face),
        zeros(n_face),
        zeros(n_face),
        zeros(n_depth, n_face),
        zeros(SVector{3, Float64}, n_face),
        zero(MVector{3, Float64}),
        zero(MVector{3, Float64}),
    )
end


# Build HierarchicalSingleAsteroidThermoPhysicalState for a HierarchicalShapeModel problem.
# Each global face with a roughness model gets its own independent sub-state.
function _build_single_state(
    problem   ::SingleAsteroidThermoPhysicalProblem{<:HierarchicalShapeModel},
    algorithm ::AbstractThermoPhysicalAlgorithm,
)
    shape   = problem.shape
    cache   = _build_cache(algorithm, problem.grid_params)
    n_depth = problem.grid_params.n_depth
    n_face  = length(shape.global_shape.faces)

    # Map each global face to an independent roughness state index (0 = no roughness)
    face_roughness_indices = zeros(Int, n_face)
    roughness_state_count  = 0
    for i in 1:n_face
        if has_roughness_model(shape, i)
            roughness_state_count += 1
            face_roughness_indices[i] = roughness_state_count
        end
    end

    # Build an independent SingleAsteroidThermoPhysicalState for each face with roughness
    roughness_states = map(findall(!=(0), face_roughness_indices)) do i
        roughness_shape = get_roughness_model(shape, i)::ShapeModel
        tp_sub = ThermoParams(
            [problem.thermo_params.conductivity[i]],
            [problem.thermo_params.density[i]],
            [problem.thermo_params.heat_capacity[i]],
            [problem.thermo_params.reflectance_vis[i]],
            [problem.thermo_params.reflectance_ir[i]],
            [problem.thermo_params.emissivity[i]],
        )
        mini_prob = SingleAsteroidThermoPhysicalProblem(
            roughness_shape, tp_sub, problem.grid_params;
            with_self_shadowing = false,
            with_self_heating   = false,
            upper_boundary_condition = problem.upper_boundary_condition,
            lower_boundary_condition = problem.lower_boundary_condition,
        )
        _build_single_state(mini_prob, algorithm)
    end

    HierarchicalSingleAsteroidThermoPhysicalState(
        problem,
        cache,
        zeros(Bool, n_face),
        zeros(n_face),
        zeros(n_face),
        zeros(n_face),
        zeros(n_depth, n_face),
        zeros(SVector{3, Float64}, n_face),
        zero(MVector{3, Float64}),
        zero(MVector{3, Float64}),
        face_roughness_indices,
        roughness_states,
    )
end


# Internal helper: build BinaryAsteroidThermoPhysicalState from a problem + algorithm.
# Both bodies always share the same algorithm and time step.
function _build_binary_state(
    problem   ::BinaryAsteroidThermoPhysicalProblem,
    algorithm ::AbstractThermoPhysicalAlgorithm,
)
    state1 = _build_single_state(problem.primary,   algorithm)
    state2 = _build_single_state(problem.secondary, algorithm)
    BinaryAsteroidThermoPhysicalState(problem, state1, state2)
end


# ╔═══════════════════════════════════════════════════════════════════╗
# ║                   Public solve interface                          ║
# ╚═══════════════════════════════════════════════════════════════════╝

"""
    solve(problem, algorithm; ephem, output, initial_temperature, show_progress=true) -> solution

Run a thermophysical simulation for a single asteroid.

# Arguments
- `problem`   : Problem definition (`SingleAsteroidThermoPhysicalProblem`)
- `algorithm` : Numerical method (`ExplicitEuler()`, `ImplicitEuler()`, or `CrankNicolson()`)

# Keyword Arguments
- `ephem`                : Ephemerides (`AbstractSingleAsteroidEphemerides`)
- `output`               : Output specification (`SingleAsteroidOutputSpec`); controls which timesteps, face indices, and physical quantities (temperatures, forces, torques) to record
- `initial_temperature`  : Initial temperature; `Real` for uniform, or `AbstractMatrix` of size `(n_depth, n_face)` [K]
- `show_progress = true` : Display progress meter during simulation

# Returns
- `SingleAsteroidThermoPhysicalSolution`

# Example
```julia
problem = SingleAsteroidThermoPhysicalProblem(shape, thermo_params;
    with_self_shadowing = true,
    with_self_heating   = true,
)
output = SingleAsteroidOutputSpec(output_times, subsurface_face_ids)
solution = solve(problem, CrankNicolson();
    ephem               = ephem,
    output              = output,
    initial_temperature = 200.0,
)
```
"""
function solve(
    problem   ::SingleAsteroidThermoPhysicalProblem,
    algorithm ::AbstractThermoPhysicalAlgorithm;
    ephem                ::AbstractSingleAsteroidEphemerides,
    output               ::SingleAsteroidOutputSpec,
    initial_temperature,
    show_progress        ::Bool = true,
)
    _validate_output_spec(output, ephem)
    _solve(problem, algorithm, ephem; output, initial_temperature, show_progress)
end


"""
    solve(problem, algorithm; ephem, output, initial_temperature_primary, initial_temperature_secondary, show_progress=true) -> solution

Run a thermophysical simulation for a binary asteroid system.

# Arguments
- `problem`   : Problem definition (`BinaryAsteroidThermoPhysicalProblem`)
- `algorithm` : Numerical method (`ExplicitEuler()`, `ImplicitEuler()`, or `CrankNicolson()`)

# Keyword Arguments
- `ephem`                         : Ephemerides (`AbstractBinaryAsteroidEphemerides`)
- `output`                        : Output specification (`BinaryAsteroidOutputSpec`); wraps a `SingleAsteroidOutputSpec` for each body
- `initial_temperature_primary`   : Initial temperature for the primary; `Real` or `AbstractMatrix` of size `(n_depth, n_face)` [K]
- `initial_temperature_secondary` : Initial temperature for the secondary; `Real` or `AbstractMatrix` of size `(n_depth, n_face)` [K]
- `show_progress = true`          : Display progress meter during simulation

# Returns
- `BinaryAsteroidThermoPhysicalSolution`

# Example
```julia
T_init = subsolar_temperature(ephem.r_sun[begin], R_vis, ε)
output = BinaryAsteroidOutputSpec(
    SingleAsteroidOutputSpec(output_times, subsurface_face_ids_pri),
    SingleAsteroidOutputSpec(output_times, subsurface_face_ids_sec),
)
solution = solve(problem, CrankNicolson();
    ephem                         = ephem,
    output                        = output,
    initial_temperature_primary   = T_init,
    initial_temperature_secondary = T_init,
)
```
"""
function solve(
    problem   ::BinaryAsteroidThermoPhysicalProblem,
    algorithm ::AbstractThermoPhysicalAlgorithm;
    ephem                        ::AbstractBinaryAsteroidEphemerides,
    output                       ::BinaryAsteroidOutputSpec,
    initial_temperature_primary,
    initial_temperature_secondary,
    show_progress                ::Bool = true,
)
    _validate_output_spec(output, ephem)
    _solve(problem, algorithm, ephem; output, initial_temperature_primary, initial_temperature_secondary, show_progress)
end


# ╔═══════════════════════════════════════════════════════════════════╗
# ║                   Internal _solve dispatch                        ║
# ╚═══════════════════════════════════════════════════════════════════╝

function _solve(
    problem   ::SingleAsteroidThermoPhysicalProblem,
    algorithm ::AbstractThermoPhysicalAlgorithm,
    ephem     ::SingleAsteroidEphemerides{Nothing};
    output              ::SingleAsteroidOutputSpec,
    initial_temperature,
    show_progress       ::Bool,
)
    state = _build_single_state(problem, algorithm)
    init_temperature!(state, initial_temperature)

    solution = SingleAsteroidThermoPhysicalSolution(state, ephem, output)

    if show_progress
        p = Progress(length(ephem.times); dt=1, desc="Running TPM...", showspeed=true)
        ProgressMeter.ijulia_behavior(:clear)
    end

    for i_time in eachindex(ephem.times)
        r☉ = ephem.r_sun[i_time]

        update_flux_all!(state, r☉)
        output.save_face_forces && update_thermal_force!(state)
        record_timestep!(solution, state, i_time)

        if show_progress
            showvalues = [
                ("Timestep     ", i_time),
                ("emitted / absorbed ", solution.emitted_power[i_time] / solution.absorbed_power[i_time]),
            ]
            ProgressMeter.next!(p; showvalues)
        end

        i_time == length(ephem.times) && break
        Δt = ephem.times[i_time+1] - ephem.times[i_time]
        update_temperature!(state, Δt)
    end

    return solution
end


function _solve(
    problem   ::SingleAsteroidThermoPhysicalProblem,
    algorithm ::AbstractThermoPhysicalAlgorithm,
    ephem     ::SingleAsteroidEphemerides{<:AbstractVector};
    output              ::SingleAsteroidOutputSpec,
    initial_temperature,
    show_progress       ::Bool,
)
    state = _build_single_state(problem, algorithm)
    init_temperature!(state, initial_temperature)

    solution = SingleAsteroidThermoPhysicalSolution(state, ephem, output)

    if show_progress
        p = Progress(length(ephem.times); dt=1, desc="Running TPM...", showspeed=true)
        ProgressMeter.ijulia_behavior(:clear)
    end

    for i_time in eachindex(ephem.times)
        r☉ = ephem.r_sun[i_time]
        R   = ephem.R_body_to_inertial[i_time]

        update_flux_all!(state, r☉)
        update_thermal_force!(state)
        record_timestep!(solution, state, i_time, R)

        if show_progress
            showvalues = [
                ("Timestep     ", i_time),
                ("emitted / absorbed ", solution.emitted_power[i_time] / solution.absorbed_power[i_time]),
            ]
            ProgressMeter.next!(p; showvalues)
        end

        i_time == length(ephem.times) && break
        Δt = ephem.times[i_time+1] - ephem.times[i_time]
        update_temperature!(state, Δt)
    end

    return solution
end


function _solve(
    problem   ::BinaryAsteroidThermoPhysicalProblem,
    algorithm ::AbstractThermoPhysicalAlgorithm,
    ephem     ::BinaryAsteroidEphemerides{Nothing};
    output                       ::BinaryAsteroidOutputSpec,
    initial_temperature_primary,
    initial_temperature_secondary,
    show_progress                ::Bool,
)
    state = _build_binary_state(problem, algorithm)
    init_temperature!(state, initial_temperature_primary, initial_temperature_secondary)

    solution = BinaryAsteroidThermoPhysicalSolution(state, ephem, output)

    if show_progress
        p = Progress(length(ephem.times); dt=1, desc="Running TPM...", showspeed=true)
        ProgressMeter.ijulia_behavior(:clear)
    end

    for i_time in eachindex(ephem.times)
        r☉₁ = ephem.r_sun[i_time]
        r₁₂ = ephem.r_secondary[i_time]
        R₁₂ = ephem.R_primary_to_secondary[i_time]

        update_flux_all!(state, r☉₁, r₁₂, R₁₂)
        (output.primary.save_face_forces || output.secondary.save_face_forces) && update_thermal_force!(state)
        record_timestep!(solution, state, i_time)

        if show_progress
            showvalues = [
                ("Timestep                   ", i_time),
                ("emitted / absorbed (primary)   ", solution.primary.emitted_power[i_time]   / solution.primary.absorbed_power[i_time]),
                ("emitted / absorbed (secondary) ", solution.secondary.emitted_power[i_time] / solution.secondary.absorbed_power[i_time]),
            ]
            ProgressMeter.next!(p; showvalues)
        end

        i_time == length(ephem.times) && break
        Δt = ephem.times[i_time+1] - ephem.times[i_time]
        update_temperature!(state, Δt)
    end

    return solution
end


function _solve(
    problem   ::BinaryAsteroidThermoPhysicalProblem,
    algorithm ::AbstractThermoPhysicalAlgorithm,
    ephem     ::BinaryAsteroidEphemerides{<:AbstractVector};
    output                       ::BinaryAsteroidOutputSpec,
    initial_temperature_primary,
    initial_temperature_secondary,
    show_progress                ::Bool,
)
    state = _build_binary_state(problem, algorithm)
    init_temperature!(state, initial_temperature_primary, initial_temperature_secondary)

    solution = BinaryAsteroidThermoPhysicalSolution(state, ephem, output)

    if show_progress
        p = Progress(length(ephem.times); dt=1, desc="Running TPM...", showspeed=true)
        ProgressMeter.ijulia_behavior(:clear)
    end

    for i_time in eachindex(ephem.times)
        r☉₁ = ephem.r_sun[i_time]
        r₁₂ = ephem.r_secondary[i_time]
        R₁₂ = ephem.R_primary_to_secondary[i_time]
        R₁ᵢ = ephem.R_primary_to_inertial[i_time]
        R₂ᵢ = R₁ᵢ * R₁₂'  # R_secondary_to_inertial = R_primary_to_inertial * R_primary_to_secondary'

        update_flux_all!(state, r☉₁, r₁₂, R₁₂)
        update_thermal_force!(state)
        record_timestep!(solution, state, i_time, R₁ᵢ, R₂ᵢ)

        if show_progress
            showvalues = [
                ("Timestep                   ", i_time),
                ("emitted / absorbed (primary)   ", solution.primary.emitted_power[i_time]   / solution.primary.absorbed_power[i_time]),
                ("emitted / absorbed (secondary) ", solution.secondary.emitted_power[i_time] / solution.secondary.absorbed_power[i_time]),
            ]
            ProgressMeter.next!(p; showvalues)
        end

        i_time == length(ephem.times) && break
        Δt = ephem.times[i_time+1] - ephem.times[i_time]
        update_temperature!(state, Δt)
    end

    return solution
end
