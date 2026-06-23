#=
tpm_solution.jl

Solution types and export functions for thermophysical simulations.

Type parameter F mirrors the ephem type parameter R:
    F = Nothing                        : forces/torques not computed (temperature only)
    F = Vector{SVector{3, Float64}}    : forces/torques stored in the inertial frame
=#

# ╔═══════════════════════════════════════════════════════════════════╗
# ║                     OutputSpec types                              ║
# ╚═══════════════════════════════════════════════════════════════════╝

"""
    struct SingleAsteroidOutputSpec

Output specification for a single-asteroid thermophysical simulation.

Encapsulates which timesteps, face indices, and physical quantities to record,
replacing the individual `times_to_save` and `face_ID` keyword arguments of `solve`.

# Fields
- `output_times`                : Timesteps at which data are saved [s]; must be a subset of `ephem.times`
- `subsurface_face_ids`         : Face indices for which to save subsurface temperature profiles
- `save_surface_temperature`    : Save surface temperature at `output_times` (default: `true`)
- `save_subsurface_temperature` : Save subsurface temperature profiles at `output_times` (default: `true`)
- `save_face_forces`            : Save per-face thermal forces at `output_times` (default: `false`)
- `save_forces`                 : Save net thermal force at `output_times` (default: `true`)
- `save_torques`                : Save net thermal torque at `output_times` (default: `true`)

# Notes
- `save_subsurface_temperature = true` requires a non-empty `subsurface_face_ids`.
- `save_forces` and `save_torques` have effect only when the ephemerides include
  `R_body_to_inertial` (i.e., `SingleAsteroidEphemerides{<:AbstractVector}`).
"""
struct SingleAsteroidOutputSpec
    output_times                ::Vector{Float64}
    subsurface_face_ids         ::Vector{Int}
    save_surface_temperature    ::Bool
    save_subsurface_temperature ::Bool
    save_face_forces            ::Bool
    save_forces                 ::Bool
    save_torques                ::Bool

    function SingleAsteroidOutputSpec(
        output_times,
        subsurface_face_ids,
        save_surface_temperature,
        save_subsurface_temperature,
        save_face_forces,
        save_forces,
        save_torques,
    )
        if save_subsurface_temperature && isempty(subsurface_face_ids)
            throw(ArgumentError(
                "subsurface_face_ids must be non-empty when save_subsurface_temperature = true"
            ))
        end
        new(
            output_times, subsurface_face_ids,
            save_surface_temperature,
            save_subsurface_temperature,
            save_face_forces,
            save_forces,
            save_torques,
        )
    end
end

SingleAsteroidOutputSpec(output_times::Vector{Float64}, subsurface_face_ids::Vector{Int}) =
    SingleAsteroidOutputSpec(output_times, subsurface_face_ids, true, true, false, true, true)

"""
    struct BinaryAsteroidOutputSpec

Output specification for a binary-asteroid thermophysical simulation.
Wraps two `SingleAsteroidOutputSpec` instances, one for each body.

# Fields
- `primary`   : Output spec for the primary body
- `secondary` : Output spec for the secondary body
"""
struct BinaryAsteroidOutputSpec
    primary   ::SingleAsteroidOutputSpec
    secondary ::SingleAsteroidOutputSpec
end

# Validate that all output_times are present in ephem.times.
function _validate_output_spec(output::SingleAsteroidOutputSpec, ephem)
    for t in output.output_times
        findfirst(isequal(t), ephem.times) === nothing &&
            throw(ArgumentError("output_times contains $t which is not in ephem.times"))
    end
end

function _validate_output_spec(output::BinaryAsteroidOutputSpec, ephem)
    _validate_output_spec(output.primary,   ephem)
    _validate_output_spec(output.secondary, ephem)
end


# ╔═══════════════════════════════════════════════════════════════════╗
# ║                     Solution types                                ║
# ╚═══════════════════════════════════════════════════════════════════╝

"""
    struct SingleAsteroidThermoPhysicalSolution{F}

Solution data for a single asteroid thermophysical simulation.

The type parameter `F` mirrors `SingleAsteroidEphemerides{R}`:
- `F = Nothing`                     : force/torque not computed; `forces = torques = nothing`
- `F = Vector{SVector{3, Float64}}` : force/torque stored in the inertial frame

# Fields
## Saved at all timesteps
- `times`           : All simulation timesteps [s]
- `absorbed_power`  : Total absorbed power on the whole surface [W]
- `emitted_power`   : Total emitted thermal radiation power from the whole surface [W]

## Saved only at `output.output_times`
- `output`                 : Output specification (controls which data are saved and when)
- `forces`                 : Net thermal force in the inertial frame [N], or `nothing`; saved when `output.save_forces = true`
- `torques`                : Net thermal torque in the inertial frame [N⋅m], or `nothing`; saved when `output.save_torques = true`
- `depth_nodes`            : Depth of each calculation node [m], size `(n_depth,)`
- `surface_temperature`    : Surface temperature [K], size `(n_face, n_save)`; saved when `output.save_surface_temperature = true`
- `subsurface_temperature` : Subsurface temperature [K] by face ID, each entry `(n_depth, n_save)`; saved when `output.save_subsurface_temperature = true`
- `face_forces`            : Per-face thermal force in the body-fixed frame [N], size `(n_face, n_save)`; saved when `output.save_face_forces = true`
"""
struct SingleAsteroidThermoPhysicalSolution{F <: Union{Nothing, AbstractVector}}
    times          ::Vector{Float64}
    absorbed_power ::Vector{Float64}
    emitted_power  ::Vector{Float64}

    output                 ::SingleAsteroidOutputSpec
    forces                 ::F
    torques                ::F
    depth_nodes            ::Vector{Float64}
    surface_temperature    ::Matrix{Float64}
    subsurface_temperature ::Dict{Int, Matrix{Float64}}
    face_forces            ::Matrix{SVector{3, Float64}}
end


"""
    struct BinaryAsteroidThermoPhysicalSolution{F}

Solution data for a binary asteroid thermophysical simulation.

# Fields
- `primary`   : Solution for the primary body
- `secondary` : Solution for the secondary body
"""
struct BinaryAsteroidThermoPhysicalSolution{F <: Union{Nothing, AbstractVector}}
    primary   ::SingleAsteroidThermoPhysicalSolution{F}
    secondary ::SingleAsteroidThermoPhysicalSolution{F}
end


# ╔═══════════════════════════════════════════════════════════════════╗
# ║                     Allocation helpers                            ║
# ╚═══════════════════════════════════════════════════════════════════╝

# Allocate a SingleAsteroidThermoPhysicalSolution{Nothing} (temperature only).
function _alloc_solution_no_force(
    state  ::SingleAsteroidThermoPhysicalState,
    times  ::Vector{Float64},
    output ::SingleAsteroidOutputSpec,
)
    n_step  = length(times)
    n_save  = length(output.output_times)
    n_face  = length(state.problem.shape.faces)
    n_depth = state.problem.thermo_params.n_depth

    absorbed_power = zeros(n_step)
    emitted_power  = zeros(n_step)
    depth_nodes            = state.problem.thermo_params.Δz * collect(0:n_depth-1)
    surface_temperature    = output.save_surface_temperature    ? zeros(n_face, n_save)                      : zeros(0, 0)
    subsurface_temperature = output.save_subsurface_temperature ? Dict{Int,Matrix{Float64}}(i => zeros(n_depth, n_save) for i in output.subsurface_face_ids) : Dict{Int,Matrix{Float64}}()
    face_forces            = output.save_face_forces            ? zeros(SVector{3,Float64}, n_face, n_save)  : zeros(SVector{3,Float64}, 0, 0)

    SingleAsteroidThermoPhysicalSolution{Nothing}(
        times, absorbed_power, emitted_power,
        output, nothing, nothing,
        depth_nodes, surface_temperature, subsurface_temperature, face_forces,
    )
end

# Allocate a SingleAsteroidThermoPhysicalSolution{Vector{SVector{3,Float64}}} (with force/torque).
function _alloc_solution_with_force(
    state  ::SingleAsteroidThermoPhysicalState,
    times  ::Vector{Float64},
    output ::SingleAsteroidOutputSpec,
)
    n_step  = length(times)
    n_save  = length(output.output_times)
    n_face  = length(state.problem.shape.faces)
    n_depth = state.problem.thermo_params.n_depth

    absorbed_power = zeros(n_step)
    emitted_power  = zeros(n_step)
    forces         = output.save_forces  ? zeros(SVector{3,Float64}, n_save) : zeros(SVector{3,Float64}, 0)
    torques        = output.save_torques ? zeros(SVector{3,Float64}, n_save) : zeros(SVector{3,Float64}, 0)
    depth_nodes            = state.problem.thermo_params.Δz * collect(0:n_depth-1)
    surface_temperature    = output.save_surface_temperature    ? zeros(n_face, n_save)                      : zeros(0, 0)
    subsurface_temperature = output.save_subsurface_temperature ? Dict{Int,Matrix{Float64}}(i => zeros(n_depth, n_save) for i in output.subsurface_face_ids) : Dict{Int,Matrix{Float64}}()
    face_forces            = output.save_face_forces            ? zeros(SVector{3,Float64}, n_face, n_save)  : zeros(SVector{3,Float64}, 0, 0)

    SingleAsteroidThermoPhysicalSolution{Vector{SVector{3,Float64}}}(
        times, absorbed_power, emitted_power,
        output, forces, torques,
        depth_nodes, surface_temperature, subsurface_temperature, face_forces,
    )
end


# ╔═══════════════════════════════════════════════════════════════════╗
# ║                     Outer constructors                            ║
# ╚═══════════════════════════════════════════════════════════════════╝

"""
    SingleAsteroidThermoPhysicalSolution(state, ephem, output)

Allocate a solution matching the type parameter of `ephem`.
"""
SingleAsteroidThermoPhysicalSolution(
    state  ::SingleAsteroidThermoPhysicalState,
    ephem  ::SingleAsteroidEphemerides{Nothing},
    output ::SingleAsteroidOutputSpec,
) = _alloc_solution_no_force(state, ephem.times, output)

SingleAsteroidThermoPhysicalSolution(
    state  ::SingleAsteroidThermoPhysicalState,
    ephem  ::SingleAsteroidEphemerides{<:AbstractVector},
    output ::SingleAsteroidOutputSpec,
) = _alloc_solution_with_force(state, ephem.times, output)

"""
    BinaryAsteroidThermoPhysicalSolution(state, ephem, output)

Allocate a binary solution matching the type parameter of `ephem`.
"""
BinaryAsteroidThermoPhysicalSolution(
    state  ::BinaryAsteroidThermoPhysicalState,
    ephem  ::BinaryAsteroidEphemerides{Nothing},
    output ::BinaryAsteroidOutputSpec,
) = BinaryAsteroidThermoPhysicalSolution{Nothing}(
    _alloc_solution_no_force(state.primary,   ephem.times, output.primary),
    _alloc_solution_no_force(state.secondary, ephem.times, output.secondary),
)

BinaryAsteroidThermoPhysicalSolution(
    state  ::BinaryAsteroidThermoPhysicalState,
    ephem  ::BinaryAsteroidEphemerides{<:AbstractVector},
    output ::BinaryAsteroidOutputSpec,
) = BinaryAsteroidThermoPhysicalSolution{Vector{SVector{3,Float64}}}(
    _alloc_solution_with_force(state.primary,   ephem.times, output.primary),
    _alloc_solution_with_force(state.secondary, ephem.times, output.secondary),
)


# ╔═══════════════════════════════════════════════════════════════════╗
# ║                     record_timestep!                              ║
# ╚═══════════════════════════════════════════════════════════════════╝

# Shared helper: record temperature and face-force data at output_times, respecting Bool flags.
function _record_temperature!(solution::SingleAsteroidThermoPhysicalSolution, state::SingleAsteroidThermoPhysicalState, i_time::Integer)
    t      = solution.times[i_time]
    i_save = something(findfirst(isequal(t), solution.output.output_times), 0)
    i_save == 0 && return

    if solution.output.save_surface_temperature
        solution.surface_temperature[:, i_save] .= surface_temperature(state)
    end

    if solution.output.save_subsurface_temperature
        for (i, T_sub) in solution.subsurface_temperature
            T_sub[:, i_save] .= state.temperature[:, i]
        end
    end

    if solution.output.save_face_forces
        solution.face_forces[:, i_save] .= state.face_forces
    end
end

"""
    record_timestep!(solution::SingleAsteroidThermoPhysicalSolution{Nothing}, state, i_time)

Record energy balance and temperature data (no force/torque).
"""
function record_timestep!(
    solution ::SingleAsteroidThermoPhysicalSolution{Nothing},
    state    ::SingleAsteroidThermoPhysicalState,
    i_time   ::Integer,
)
    solution.absorbed_power[i_time] = energy_in(state)
    solution.emitted_power[i_time]  = energy_out(state)
    _record_temperature!(solution, state, i_time)
end

"""
    record_timestep!(solution::SingleAsteroidThermoPhysicalSolution{<:AbstractVector}, state, i_time, R)

Record energy balance at all timesteps, and at `output_times`: force/torque rotated to the
inertial frame via `R` and temperature data.
"""
function record_timestep!(
    solution ::SingleAsteroidThermoPhysicalSolution{<:AbstractVector},
    state    ::SingleAsteroidThermoPhysicalState,
    i_time   ::Integer,
    R        ::SMatrix{3,3,Float64,9},
)
    solution.absorbed_power[i_time] = energy_in(state)
    solution.emitted_power[i_time]  = energy_out(state)

    i_save = something(findfirst(isequal(solution.times[i_time]), solution.output.output_times), 0)
    if i_save != 0
        solution.output.save_forces  && (solution.forces[i_save]  = R * state.force)
        solution.output.save_torques && (solution.torques[i_save] = R * state.torque)
    end
    _record_temperature!(solution, state, i_time)
end

"""
    record_timestep!(solution::BinaryAsteroidThermoPhysicalSolution{Nothing}, state, i_time)
"""
function record_timestep!(
    solution ::BinaryAsteroidThermoPhysicalSolution{Nothing},
    state    ::BinaryAsteroidThermoPhysicalState,
    i_time   ::Integer,
)
    record_timestep!(solution.primary,   state.primary,   i_time)
    record_timestep!(solution.secondary, state.secondary, i_time)
end

"""
    record_timestep!(solution::BinaryAsteroidThermoPhysicalSolution{<:AbstractVector}, state, i_time, R_p2i, R_s2i)
"""
function record_timestep!(
    solution ::BinaryAsteroidThermoPhysicalSolution{<:AbstractVector},
    state    ::BinaryAsteroidThermoPhysicalState,
    i_time   ::Integer,
    R₁ᵢ      ::SMatrix{3,3,Float64,9},
    R₂ᵢ      ::SMatrix{3,3,Float64,9},
)
    record_timestep!(solution.primary,   state.primary,   i_time, R₁ᵢ)
    record_timestep!(solution.secondary, state.secondary, i_time, R₂ᵢ)
end


# ╔═══════════════════════════════════════════════════════════════════╗
# ║                     export_solution                               ║
# ╚═══════════════════════════════════════════════════════════════════╝

function _export_diagnostics(dirpath, solution::SingleAsteroidThermoPhysicalSolution)
    df = DataFrame(
        time           = solution.times,
        absorbed_power = solution.absorbed_power,
        emitted_power  = solution.emitted_power,
    )

    filepath = joinpath(dirpath, "diagnostics.csv")
    CSV.write(filepath, df)
end

function _export_thermal_net_forces(dirpath, solution::SingleAsteroidThermoPhysicalSolution{<:AbstractVector})
    df = DataFrame(time = solution.output.output_times)
    if solution.output.save_forces
        df.force_x  = [f[1] for f in solution.forces]
        df.force_y  = [f[2] for f in solution.forces]
        df.force_z  = [f[3] for f in solution.forces]
    end
    if solution.output.save_torques
        df.torque_x = [τ[1] for τ in solution.torques]
        df.torque_y = [τ[2] for τ in solution.torques]
        df.torque_z = [τ[3] for τ in solution.torques]
    end

    filepath = joinpath(dirpath, "thermal_net_forces.csv")
    CSV.write(filepath, df)
end

function _export_surface_temperature(dirpath, solution::SingleAsteroidThermoPhysicalSolution)
    df = hcat(
        DataFrame(time = solution.output.output_times),
        DataFrame(
            solution.surface_temperature',
            ["face_$(i)" for i in 1:size(solution.surface_temperature, 1)],
        ),
    )

    filepath = joinpath(dirpath, "surface_temperature.csv")
    CSV.write(filepath, df)
end

function _export_subsurface_temperature(dirpath, solution::SingleAsteroidThermoPhysicalSolution)
    output_times = solution.output.output_times
    nrows = length(solution.depth_nodes) * length(output_times)
    df = DataFrame(
        time  = reshape([t for _ in solution.depth_nodes, t in output_times], nrows),
        depth = reshape([d for d in solution.depth_nodes, _ in output_times], nrows),
    )
    for (i, T_sub) in collect(solution.subsurface_temperature)
        df[:, "face_$(i)"] = reshape(T_sub, length(T_sub))
    end
    keys_sorted = sort(names(df[:, 3:end]), by=x->parse(Int, replace(x, r"[^0-9]" => "")))
    df = df[:, ["time", "depth", keys_sorted...]]
    CSV.write(joinpath(dirpath, "subsurface_temperature.csv"), df)
end

function _export_thermal_face_forces(dirpath, solution::SingleAsteroidThermoPhysicalSolution)
    output_times = solution.output.output_times
    n_face = size(solution.face_forces, 1)
    n_step = size(solution.face_forces, 2)
    nrows  = n_face * n_step
    df = DataFrame(
        time    = reshape([t for _ in 1:n_face, t in output_times], nrows),
        face_id = reshape([i for i in 1:n_face, _ in output_times], nrows),
    )
    df.force_x = reshape([f[1] for f in solution.face_forces], nrows)
    df.force_y = reshape([f[2] for f in solution.face_forces], nrows)
    df.force_z = reshape([f[3] for f in solution.face_forces], nrows)

    filepath = joinpath(dirpath, "thermal_face_forces.csv") 
    CSV.write(filepath, df)
end


"""
    export_solution(dirpath, solution::SingleAsteroidThermoPhysicalSolution{Nothing})

Export simulation results to CSV files in `dirpath`.

Files written depend on the `output` specification:
- `diagnostics.csv`            : always (`absorbed_power`, `emitted_power` at all timesteps)
- `surface_temperature.csv`    : when `output.save_surface_temperature = true`
- `subsurface_temperature.csv` : when `output.save_subsurface_temperature = true`
- `thermal_face_forces.csv`    : when `output.save_face_forces = true`
"""
function export_solution(dirpath, solution::SingleAsteroidThermoPhysicalSolution{Nothing})
    _export_diagnostics(dirpath, solution)
    solution.output.save_surface_temperature    && _export_surface_temperature(dirpath, solution)
    solution.output.save_subsurface_temperature && _export_subsurface_temperature(dirpath, solution)
    solution.output.save_face_forces            && _export_thermal_face_forces(dirpath, solution)
end


"""
    export_solution(dirpath, solution::SingleAsteroidThermoPhysicalSolution{<:AbstractVector})

Export simulation results to CSV files in `dirpath`.

Files written depend on the `output` specification:
- `diagnostics.csv`            : always (`absorbed_power`, `emitted_power` at all timesteps)
- `thermal_net_forces.csv`     : when `output.save_forces = true` or `output.save_torques = true`
- `surface_temperature.csv`    : when `output.save_surface_temperature = true`
- `subsurface_temperature.csv` : when `output.save_subsurface_temperature = true`
- `thermal_face_forces.csv`    : when `output.save_face_forces = true`
"""
function export_solution(dirpath, solution::SingleAsteroidThermoPhysicalSolution{<:AbstractVector})
    _export_diagnostics(dirpath, solution)
    solution.output.save_surface_temperature    && _export_surface_temperature(dirpath, solution)
    solution.output.save_subsurface_temperature && _export_subsurface_temperature(dirpath, solution)
    solution.output.save_face_forces            && _export_thermal_face_forces(dirpath, solution)

    if solution.output.save_forces || solution.output.save_torques
        _export_thermal_net_forces(dirpath, solution)
    end
end


"""
    export_solution(dirpath, solution::BinaryAsteroidThermoPhysicalSolution)

Export results for both bodies to `dirpath/primary/` and `dirpath/secondary/`.
"""
function export_solution(dirpath, solution::BinaryAsteroidThermoPhysicalSolution)
    dirpath1 = joinpath(dirpath, "primary")
    dirpath2 = joinpath(dirpath, "secondary")
    mkpath(dirpath1)
    mkpath(dirpath2)
    export_solution(dirpath1, solution.primary)
    export_solution(dirpath2, solution.secondary)
end
