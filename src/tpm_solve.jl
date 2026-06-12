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
    solve(problem, algorithm; ephem, times_to_save, face_ID, T₀, show_progress=true) -> solution

Run a thermophysical simulation for a single asteroid.

# Arguments
- `problem`   : Problem definition (`SingleAsteroidThermoPhysicalProblem`)
- `algorithm` : Numerical method (`ExplicitEuler()`, `ImplicitEuler()`, or `CrankNicolson()`)

# Keyword Arguments
- `ephem`                : Ephemerides (`AbstractSingleAsteroidEphemerides`)
- `times_to_save`        : Time points at which to save detailed temperature data [s]
- `face_ID`              : Face indices for which to save subsurface temperature profiles
- `T₀`                   : Initial temperature; `Real` for uniform, or `AbstractMatrix` of size `(n_depth, n_face)` [K]
- `show_progress = true` : Display progress meter during simulation

# Returns
- `SingleAsteroidThermoPhysicalSolution`

# Example
```julia
problem = SingleAsteroidThermoPhysicalProblem(shape, thermo_params;
    with_self_shadowing = true,
    with_self_heating   = true,
)
solution = solve(problem, CrankNicolson();
    ephem         = ephem,
    times_to_save = times_to_save,
    face_ID       = face_ID,
    T₀            = 200.0,
)
```
"""
function solve(
    problem   ::SingleAsteroidThermoPhysicalProblem,
    algorithm ::AbstractThermoPhysicalAlgorithm;
    ephem         ::AbstractSingleAsteroidEphemerides,
    times_to_save ::Vector{Float64},
    face_ID       ::Vector{Int},
    T₀,
    show_progress ::Bool = true,
)
    _solve(problem, algorithm, ephem; times_to_save, face_ID, T₀, show_progress)
end


"""
    solve(problem, algorithm; ephem, times_to_save, face_ID_pri, face_ID_sec, T₀_primary, T₀_secondary, show_progress=true) -> solution

Run a thermophysical simulation for a binary asteroid system.

# Arguments
- `problem`   : Problem definition (`BinaryAsteroidThermoPhysicalProblem`)
- `algorithm` : Numerical method (`ExplicitEuler()`, `ImplicitEuler()`, or `CrankNicolson()`)

# Keyword Arguments
- `ephem`                  : Ephemerides (`AbstractBinaryAsteroidEphemerides`)
- `times_to_save`          : Time points at which to save detailed temperature data [s]
- `face_ID_pri`            : Face indices for subsurface profiles of the primary
- `face_ID_sec`            : Face indices for subsurface profiles of the secondary
- `T₀_primary`             : Initial temperature for the primary; `Real` or `AbstractMatrix` of size `(n_depth, n_face)` [K]
- `T₀_secondary`           : Initial temperature for the secondary; `Real` or `AbstractMatrix` of size `(n_depth, n_face)` [K]
- `show_progress = true`   : Display progress meter during simulation

# Returns
- `BinaryAsteroidThermoPhysicalSolution`

# Example
```julia
T₀_primary   = subsolar_temperature(ephem.r_sun[begin], params1)
T₀_secondary = subsolar_temperature(ephem.r_sun[begin], params2)
solution = solve(problem, CrankNicolson();
    ephem         = ephem,
    times_to_save = times_to_save,
    face_ID_pri   = face_ID_pri,
    face_ID_sec   = face_ID_sec,
    T₀_primary    = T₀_primary,
    T₀_secondary  = T₀_secondary,
)
```
"""
function solve(
    problem   ::BinaryAsteroidThermoPhysicalProblem,
    algorithm ::AbstractThermoPhysicalAlgorithm;
    ephem         ::AbstractBinaryAsteroidEphemerides,
    times_to_save ::Vector{Float64},
    face_ID_pri   ::Vector{Int},
    face_ID_sec   ::Vector{Int},
    T₀_primary,
    T₀_secondary,
    show_progress ::Bool = true,
)
    _solve(problem, algorithm, ephem; times_to_save, face_ID_pri, face_ID_sec, T₀_primary, T₀_secondary, show_progress)
end


# ╔═══════════════════════════════════════════════════════════════════╗
# ║                   Internal _solve dispatch                        ║
# ╚═══════════════════════════════════════════════════════════════════╝

function _solve(
    problem   ::SingleAsteroidThermoPhysicalProblem,
    algorithm ::AbstractThermoPhysicalAlgorithm,
    ephem     ::SingleAsteroidEphemerides;
    times_to_save ::Vector{Float64},
    face_ID       ::Vector{Int},
    T₀,
    show_progress ::Bool,
)
    state = _build_single_state(problem, algorithm)
    init_temperature!(state, T₀)

    solution = SingleAsteroidThermoPhysicalSolution(state, ephem, times_to_save, face_ID)

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
                ("E_out / E_in ", solution.E_out[i_time] / solution.E_in[i_time]),
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
    ephem     ::SingleAsteroidEphemeridesWithDynamics;
    kwargs...,
)
    error("solve with SingleAsteroidEphemeridesWithDynamics is not yet implemented")
end


function _solve(
    problem   ::BinaryAsteroidThermoPhysicalProblem,
    algorithm ::AbstractThermoPhysicalAlgorithm,
    ephem     ::BinaryAsteroidEphemerides;
    times_to_save ::Vector{Float64},
    face_ID_pri   ::Vector{Int},
    face_ID_sec   ::Vector{Int},
    T₀_primary,
    T₀_secondary,
    show_progress ::Bool,
)
    state = _build_binary_state(problem, algorithm)
    init_temperature!(state, T₀_primary, T₀_secondary)

    solution = BinaryAsteroidThermoPhysicalSolution(state, ephem, times_to_save, face_ID_pri, face_ID_sec)

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
                ("E_out / E_in for primary   ", solution.primary.E_out[i_time] / solution.primary.E_in[i_time]),
                ("E_out / E_in for secondary ", solution.secondary.E_out[i_time] / solution.secondary.E_in[i_time]),
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
    ephem     ::BinaryAsteroidEphemeridesWithDynamics;
    kwargs...,
)
    error("solve with BinaryAsteroidEphemeridesWithDynamics is not yet implemented")
end
