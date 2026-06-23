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
_build_cache(::ExplicitEuler, thermo_params) = ExplicitEulerCache(thermo_params)
_build_cache(::ImplicitEuler, thermo_params) = ImplicitEulerCache(thermo_params)
_build_cache(::CrankNicolson, thermo_params) = CrankNicolsonCache(thermo_params)


# Internal helper: build SingleAsteroidThermoPhysicalState from a problem + algorithm.
# broadcast_thermo_params! has already been called in the problem constructor.
function _build_single_state(
    problem   ::SingleAsteroidThermoPhysicalProblem,
    algorithm ::AbstractThermoPhysicalAlgorithm,
)
    cache   = _build_cache(algorithm, problem.thermo_params)
    n_depth = problem.thermo_params.n_depth
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
- `ephem`                   : Ephemerides (`AbstractSingleAsteroidEphemerides`)
- `output`                  : Output specification (`SingleAsteroidOutputSpec`); defines which timesteps and face indices to record in detail
- `initial_temperature`     : Initial temperature; `Real` for uniform, or `AbstractMatrix` of size `(n_depth, n_face)` [K]
- `show_progress = true`    : Display progress meter during simulation

# Returns
- `SingleAsteroidThermoPhysicalSolution`

# Example
```julia
problem = SingleAsteroidThermoPhysicalProblem(shape, thermo_params;
    with_self_shadowing = true,
    with_self_heating   = true,
)
output = SingleAsteroidOutputSpec(times_to_save, face_ID)
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
    SingleAsteroidOutputSpec(times_to_save, face_ID_pri),
    SingleAsteroidOutputSpec(times_to_save, face_ID_sec),
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
        update_thermal_force!(state)
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
        update_thermal_force!(state)
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
