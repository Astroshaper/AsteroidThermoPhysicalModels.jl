#=
tpm_solve.jl

Implementation of the solve interface for thermophysical problems.
Extends CommonSolve.solve for AbstractThermoPhysicalProblem types.
=#

# Build the appropriate HeatConductionCache from an algorithm type
_build_cache(::ExplicitEuler, thermo_params) = ExplicitEulerCache(thermo_params)
_build_cache(::ImplicitEuler, thermo_params) = ImplicitEulerCache(thermo_params)
_build_cache(::CrankNicolson, thermo_params) = CrankNicolsonCache(thermo_params)


# Internal helper: build SingleAsteroidThermoPhysicalModel from a problem + algorithm,
# bypassing broadcast_thermo_params! (already called in the problem constructor).
function _build_single_tpm(
    problem   ::SingleAsteroidThermoPhysicalProblem,
    algorithm ::AbstractThermoPhysicalAlgorithm,
)
    cache   = _build_cache(algorithm, problem.thermo_params)
    n_depth = problem.thermo_params.n_depth
    n_face  = length(problem.shape.faces)

    SingleAsteroidThermoPhysicalModel(
        problem.shape,
        problem.thermo_params,
        zeros(Bool, n_face),
        zeros(n_face),
        zeros(n_face),
        zeros(n_face),
        zeros(n_depth, n_face),
        zeros(SVector{3, Float64}, n_face),
        zero(MVector{3, Float64}),
        zero(MVector{3, Float64}),
        problem.with_self_shadowing,
        problem.with_self_heating,
        cache,
        problem.upper_boundary_condition,
        problem.lower_boundary_condition,
    )
end


"""
    solve(problem, algorithm; ephem, times_to_save, face_ID, T₀=200.0, show_progress=true) -> solution

Run a thermophysical simulation for a single asteroid.

# Arguments
- `problem`   : Problem definition (`SingleAsteroidThermoPhysicalProblem`)
- `algorithm` : Numerical method (`ExplicitEuler()`, `ImplicitEuler()`, or `CrankNicolson()`)

# Keyword Arguments
- `ephem`                : Ephemerides with fields `time::Vector{Float64}` and `sun::Vector{SVector{3}}`
- `times_to_save`        : Time points at which to save detailed temperature data [s]
- `face_ID`              : Face indices for which to save subsurface temperature profiles
- `T₀ = 200.0`           : Initial temperature for all cells [K]
- `show_progress = true` : Display progress meter during simulation

# Returns
- `SingleAsteroidThermoPhysicalModelResult`

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
    ephem,
    times_to_save ::Vector{Float64},
    face_ID       ::Vector{Int},
    T₀            ::Real = 200.0,
    show_progress ::Bool = true,
)
    stpm = _build_single_tpm(problem, algorithm)
    init_temperature!(stpm, T₀)

    solution = SingleAsteroidThermoPhysicalModelResult(stpm, ephem, times_to_save, face_ID)

    if show_progress
        p = Progress(length(ephem.time); dt=1, desc="Running TPM...", showspeed=true)
        ProgressMeter.ijulia_behavior(:clear)
    end

    for i_time in eachindex(ephem.time)
        r☉ = ephem.sun[i_time]

        update_flux_all!(stpm, r☉)
        update_thermal_force!(stpm)
        update_TPM_result!(solution, stpm, i_time)

        if show_progress
            showvalues = [
                ("Timestep     ", i_time),
                ("E_out / E_in ", solution.E_out[i_time] / solution.E_in[i_time]),
            ]
            ProgressMeter.next!(p; showvalues)
        end

        i_time == length(ephem.time) && break
        Δt = ephem.time[i_time+1] - ephem.time[i_time]
        update_temperature!(stpm, Δt)
    end

    return solution
end


"""
    solve(problem, algorithm; ephem, times_to_save, face_ID_pri, face_ID_sec, T₀=200.0, show_progress=true) -> solution

Run a thermophysical simulation for a binary asteroid system.

# Arguments
- `problem`   : Problem definition (`BinaryAsteroidThermoPhysicalProblem`)
- `algorithm` : Numerical method (`ExplicitEuler()`, `ImplicitEuler()`, or `CrankNicolson()`)

# Keyword Arguments
- `ephem`               : Ephemerides with fields `time`, `sun`, `sec`, `P2S`
- `times_to_save`       : Time points at which to save detailed temperature data [s]
- `face_ID_pri`         : Face indices for subsurface profiles of the primary
- `face_ID_sec`         : Face indices for subsurface profiles of the secondary
- `T₀ = 200.0`          : Initial temperature for all cells in both bodies [K]
- `show_progress = true` : Display progress meter during simulation

# Returns
- `BinaryAsteroidThermoPhysicalModelResult`

# Example
```julia
solution = solve(problem, CrankNicolson();
    ephem         = ephem,
    times_to_save = times_to_save,
    face_ID_pri   = face_ID_pri,
    face_ID_sec   = face_ID_sec,
    T₀            = 200.0,
)
```
"""
function solve(
    problem   ::BinaryAsteroidThermoPhysicalProblem,
    algorithm ::AbstractThermoPhysicalAlgorithm;
    ephem,
    times_to_save ::Vector{Float64},
    face_ID_pri   ::Vector{Int},
    face_ID_sec   ::Vector{Int},
    T₀            ::Real = 200.0,
    show_progress ::Bool = true,
)
    stpm1 = _build_single_tpm(problem.primary,   algorithm)
    stpm2 = _build_single_tpm(problem.secondary, algorithm)

    btpm = BinaryAsteroidThermoPhysicalModel(
        stpm1, stpm2,
        problem.with_mutual_shadowing,
        problem.with_mutual_heating,
    )
    init_temperature!(btpm, T₀)

    solution = BinaryAsteroidThermoPhysicalModelResult(btpm, ephem, times_to_save, face_ID_pri, face_ID_sec)

    if show_progress
        p = Progress(length(ephem.time); dt=1, desc="Running TPM...", showspeed=true)
        ProgressMeter.ijulia_behavior(:clear)
    end

    for i_time in eachindex(ephem.time)
        r☉₁ = ephem.sun[i_time]
        r₁₂ = ephem.sec[i_time]
        R₁₂ = ephem.P2S[i_time]

        update_flux_all!(btpm, r☉₁, r₁₂, R₁₂)
        update_thermal_force!(btpm)
        update_TPM_result!(solution, btpm, i_time)

        if show_progress
            showvalues = [
                ("Timestep                   ", i_time),
                ("E_out / E_in for primary   ", solution.pri.E_out[i_time] / solution.pri.E_in[i_time]),
                ("E_out / E_in for secondary ", solution.sec.E_out[i_time] / solution.sec.E_in[i_time]),
            ]
            ProgressMeter.next!(p; showvalues)
        end

        i_time == length(ephem.time) && break
        Δt = ephem.time[i_time+1] - ephem.time[i_time]
        update_temperature!(btpm, Δt)
    end

    return solution
end
