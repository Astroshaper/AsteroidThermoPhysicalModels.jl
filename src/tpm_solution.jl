#=
tpm_solution.jl

Solution types and export functions for thermophysical simulations.
=#

# ╔═══════════════════════════════════════════════════════════════════╗
# ║                     OutputSpec types                              ║
# ╚═══════════════════════════════════════════════════════════════════╝

"""
    struct SingleAsteroidOutputSpec

Output specification for a single-asteroid thermophysical simulation.

Encapsulates which timesteps, face indices, and physical quantities to record. Quantities
recorded for every face are switched on and off by a `save_*` flag; quantities recorded for
selected faces only are switched on by listing those faces — an empty list means "not saved".

# Fields
- `output_times`             : Timesteps at which data are saved [s]; must be a subset of `ephem.times`
- `save_surface_temperature` : Save the surface temperature of every face at `output_times` (default: `true`)
- `subsurface_face_ids`      : Faces whose subsurface temperature profiles to save at `output_times` (default: none)
- `roughness_face_ids`       : Faces whose roughness-model surface temperatures (every sub-face) to save
  at `output_times` (default: none); requires a shape with surface roughness
- `save_face_forces`         : Save the thermal force on every face at `output_times` (default: `false`)
- `save_forces`              : Save the net thermal force at `output_times` (default: `false`)
- `save_torques`             : Save the net thermal torque at `output_times` (default: `false`)

# Notes
- `roughness_face_ids` requires a problem built on a shape with surface roughness and a roughness
  model on every listed face; both are checked when the solution is allocated at `solve` time.
- `save_face_forces` stores per-face forces in the body-fixed frame; it works with both
  `SingleAsteroidEphemerides{Nothing}` and `SingleAsteroidEphemerides{<:AbstractVector}`.
- `save_forces` and `save_torques` require ephemerides with `R_body_to_inertial`
  (i.e., `SingleAsteroidEphemerides{<:AbstractVector}`); using them with rotation-free
  ephemerides raises an `ArgumentError` at `solve` time.

# Example
```julia
output = SingleAsteroidOutputSpec(output_times;
    save_surface_temperature = true,
    subsurface_face_ids      = [1, 2, 3],
    roughness_face_ids       = [1, 7],
    save_face_forces         = false,
    save_forces              = true,
    save_torques             = true,
)
```
"""
struct SingleAsteroidOutputSpec
    output_times             ::Vector{Float64}
    save_surface_temperature ::Bool
    subsurface_face_ids      ::Vector{Int}
    roughness_face_ids       ::Vector{Int}
    save_face_forces         ::Bool
    save_forces              ::Bool
    save_torques             ::Bool
end

"""
    SingleAsteroidOutputSpec(output_times; kwargs...)

Keyword constructor; see [`SingleAsteroidOutputSpec`](@ref) for the meaning of each keyword.
"""
function SingleAsteroidOutputSpec(
    output_times             ::AbstractVector{<:Real};
    save_surface_temperature ::Bool = true,
    subsurface_face_ids      ::AbstractVector{<:Integer} = Int[],
    roughness_face_ids       ::AbstractVector{<:Integer} = Int[],
    save_face_forces         ::Bool = false,
    save_forces              ::Bool = false,
    save_torques             ::Bool = false,
)
    SingleAsteroidOutputSpec(
        convert(Vector{Float64}, output_times),
        save_surface_temperature,
        convert(Vector{Int}, subsurface_face_ids),
        convert(Vector{Int}, roughness_face_ids),
        save_face_forces,
        save_forces,
        save_torques,
    )
end

"""
    struct BinaryAsteroidOutputSpec

Output specification for a binary-asteroid thermophysical simulation.
Wraps two `SingleAsteroidOutputSpec` instances, one for each body.

# Fields
- `primary`   : Output spec for the primary body
- `secondary` : Output spec for the secondary body

# Example
```julia
# Convenience constructor: shared Bool flags, separate output_times and subsurface_face_ids
output = BinaryAsteroidOutputSpec(output_times_primary, output_times_secondary;
    subsurface_face_ids_primary   = [1, 2],
    subsurface_face_ids_secondary = [3],
    save_surface_temperature = true,
    save_face_forces         = false,
    save_forces              = true,
    save_torques             = true,
)

# Base constructor: independent spec per body
output_primary   = SingleAsteroidOutputSpec(output_times_primary;   subsurface_face_ids=[1, 2], save_forces=true)
output_secondary = SingleAsteroidOutputSpec(output_times_secondary; subsurface_face_ids=[3])
output = BinaryAsteroidOutputSpec(output_primary, output_secondary)
```
"""
struct BinaryAsteroidOutputSpec
    primary   ::SingleAsteroidOutputSpec
    secondary ::SingleAsteroidOutputSpec
end

function BinaryAsteroidOutputSpec(
    output_times_primary   ::AbstractVector{<:Real},
    output_times_secondary ::AbstractVector{<:Real};
    subsurface_face_ids_primary   ::AbstractVector{<:Integer} = Int[],
    subsurface_face_ids_secondary ::AbstractVector{<:Integer} = Int[],
    save_surface_temperature ::Bool = true,
    save_face_forces         ::Bool = false,
    save_forces              ::Bool = false,
    save_torques             ::Bool = false,
)
    output_primary = SingleAsteroidOutputSpec(output_times_primary;
        subsurface_face_ids = subsurface_face_ids_primary,
        save_surface_temperature, save_face_forces, save_forces, save_torques,
    )
    output_secondary = SingleAsteroidOutputSpec(output_times_secondary;
        subsurface_face_ids = subsurface_face_ids_secondary,
        save_surface_temperature, save_face_forces, save_forces, save_torques,
    )
    BinaryAsteroidOutputSpec(output_primary, output_secondary)
end


# ╔═══════════════════════════════════════════════════════════════════╗
# ║                  OutputSpec validation                            ║
# ╚═══════════════════════════════════════════════════════════════════╝

# Validate that all output_times are present in ephem.times.
function _validate_output_spec(output::SingleAsteroidOutputSpec, ephem)
    for t in output.output_times
        findfirst(isequal(t), ephem.times) === nothing &&
            throw(ArgumentError("output_times contains $t which is not in ephem.times"))
    end
end

# For {Nothing} ephem, forces and torques cannot be computed (no rotation matrices).
function _validate_output_spec(output::SingleAsteroidOutputSpec, ephem::SingleAsteroidEphemerides{Nothing})
    for t in output.output_times
        findfirst(isequal(t), ephem.times) === nothing &&
            throw(ArgumentError("output_times contains $t which is not in ephem.times"))
    end

    if output.save_forces || output.save_torques
        throw(ArgumentError(
            "save_forces and save_torques require R_body_to_inertial in ephemerides " *
            "(use SingleAsteroidEphemerides with rotation matrices)"
        ))
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
    struct SingleAsteroidThermoPhysicalSolution

Solution data for a single asteroid thermophysical simulation.

# Fields
## Saved at all timesteps
- `times`           : All simulation timesteps [s]
- `absorbed_power`  : Total absorbed power on the whole surface [W]
- `emitted_power`   : Total emitted thermal radiation power from the whole surface [W]

## Metadata
- `output`      : Output specification (controls which data are saved and when)
- `depth_nodes` : Depth of each subsurface calculation node [m], size `(n_depth,)`

## Saved only at `output.output_times` (`nothing` when the corresponding flag is `false`)
- `surface_temperature`    : Surface temperature [K], size `(n_face, n_save)`, or `nothing`
- `subsurface_temperature` : Subsurface temperature [K] by face ID, each entry `(n_depth, n_save)`, or `nothing`
- `face_forces`            : Per-face thermal force in the body-fixed frame [N], size `(n_face, n_save)`, or `nothing`
- `forces`                 : Net thermal force in the inertial frame [N], size `(n_save,)`, or `nothing`
- `torques`                : Net thermal torque in the inertial frame [N⋅m], size `(n_save,)`, or `nothing`
- `roughness_surface_temperature` : Surface temperature of the sub-faces of the roughness model, by
  global face ID, each entry `(n_sub, n_save)` [K], or `nothing`. Sub-face `j` is face `j` of the
  roughness model attached to that global face

# Notes
- `forces` and `torques` are non-`nothing` only when the ephemerides include
  `R_body_to_inertial` (i.e., `SingleAsteroidEphemerides{<:AbstractVector}`) and the
  corresponding flag in `output` is `true`.
- `surface_temperature` and `subsurface_temperature` are those of the global faces, i.e. the
  smooth-surface baseline of a hierarchical run; `roughness_surface_temperature` is where the
  rough-surface temperatures live.
"""
struct SingleAsteroidThermoPhysicalSolution
    times          ::Vector{Float64}
    absorbed_power ::Vector{Float64}
    emitted_power  ::Vector{Float64}

    output                 ::SingleAsteroidOutputSpec
    depth_nodes            ::Vector{Float64}
    surface_temperature    ::Union{Nothing, Matrix{Float64}}
    subsurface_temperature ::Union{Nothing, Dict{Int, Matrix{Float64}}}
    face_forces            ::Union{Nothing, Matrix{SVector{3, Float64}}}
    forces                 ::Union{Nothing, Vector{SVector{3, Float64}}}
    torques                ::Union{Nothing, Vector{SVector{3, Float64}}}
    roughness_surface_temperature ::Union{Nothing, Dict{Int, Matrix{Float64}}}
end


"""
    struct BinaryAsteroidThermoPhysicalSolution

Solution data for a binary asteroid thermophysical simulation.

# Fields
- `primary`   : Solution for the primary body
- `secondary` : Solution for the secondary body
"""
struct BinaryAsteroidThermoPhysicalSolution
    primary   ::SingleAsteroidThermoPhysicalSolution
    secondary ::SingleAsteroidThermoPhysicalSolution
end


# ╔═══════════════════════════════════════════════════════════════════╗
# ║                  Solution constructors                            ║
# ╚═══════════════════════════════════════════════════════════════════╝

function _build_single_solution(
    state  ::SingleAsteroidThermoPhysicalState,
    times  ::Vector{Float64},
    output ::SingleAsteroidOutputSpec,
)
    n_step  = length(times)
    n_save  = length(output.output_times)
    n_face  = size(state.temperature, 2)
    n_depth = state.problem.grid_params.n_depth

    absorbed_power = zeros(n_step)
    emitted_power  = zeros(n_step)

    depth_nodes            = state.problem.grid_params.Δz * collect(0:n_depth-1)
    surface_temperature    = output.save_surface_temperature    ? zeros(n_face, n_save) : nothing
    subsurface_temperature = !isempty(output.subsurface_face_ids) ? Dict{Int,Matrix{Float64}}(i => zeros(n_depth, n_save) for i in output.subsurface_face_ids) : nothing
    face_forces            = output.save_face_forces            ? zeros(SVector{3,Float64}, n_face, n_save) : nothing
    forces                 = output.save_forces  ? zeros(SVector{3,Float64}, n_save) : nothing
    torques                = output.save_torques ? zeros(SVector{3,Float64}, n_save) : nothing
    roughness_surface_temperature = !isempty(output.roughness_face_ids) ?
        Dict{Int,Matrix{Float64}}(i => zeros(_n_roughness_faces(state, i), n_save) for i in output.roughness_face_ids) : nothing

    SingleAsteroidThermoPhysicalSolution(
        times, absorbed_power, emitted_power,
        output,
        depth_nodes, surface_temperature, subsurface_temperature, face_forces,
        forces, torques,
        roughness_surface_temperature,
    )
end

# Number of sub-faces of the roughness model on face `i`; raises when the state has no
# roughness models at all, or face `i` carries none, since there would be nothing to save.
function _n_roughness_faces(state::SingleAsteroidThermoPhysicalState, i::Integer)
    isempty(state.roughness_states) && throw(ArgumentError(
        "roughness_face_ids requires a shape with surface roughness (see `add_roughness_models!`)"
    ))
    k = state.face_roughness_indices[i]
    k == 0 && throw(ArgumentError(
        "roughness_face_ids contains face $i, which has no roughness model"
    ))
    return size(state.roughness_states[k].temperature, 2)
end

"""
    SingleAsteroidThermoPhysicalSolution(state, ephem, output)

Allocate a solution from the given state, ephemerides, and output specification.
"""
SingleAsteroidThermoPhysicalSolution(
    state  ::SingleAsteroidThermoPhysicalState,
    ephem  ::AbstractSingleAsteroidEphemerides,
    output ::SingleAsteroidOutputSpec,
) = _build_single_solution(state, ephem.times, output)

"""
    BinaryAsteroidThermoPhysicalSolution(state, ephem, output)

Allocate a binary solution from the given state, ephemerides, and output specification.
"""
BinaryAsteroidThermoPhysicalSolution(
    state  ::BinaryAsteroidThermoPhysicalState,
    ephem  ::AbstractBinaryAsteroidEphemerides,
    output ::BinaryAsteroidOutputSpec,
) = BinaryAsteroidThermoPhysicalSolution(
    _build_single_solution(state.primary,   ephem.times, output.primary),
    _build_single_solution(state.secondary, ephem.times, output.secondary),
)


# ╔═══════════════════════════════════════════════════════════════════╗
# ║                     record_timestep!                              ║
# ╚═══════════════════════════════════════════════════════════════════╝

"""
    record_timestep!(solution, state, i_time)

Record simulation data for a single asteroid at timestep `i_time`.

Saves `absorbed_power` and `emitted_power` at every timestep. At timesteps that coincide
with `solution.output.output_times`, also records snapshot data (surface temperature,
subsurface temperature, face forces) according to the flags in `solution.output`.

# Arguments
- `solution` : Solution container (`SingleAsteroidThermoPhysicalSolution`)
- `state`    : Current simulation state (`SingleAsteroidThermoPhysicalState`)
- `i_time`   : Index into `solution.times` for the current timestep
"""
function record_timestep!(
    solution ::SingleAsteroidThermoPhysicalSolution,
    state    ::SingleAsteroidThermoPhysicalState,
    i_time   ::Integer,
)
    solution.absorbed_power[i_time] = integrate_absorbed_power(state)
    solution.emitted_power[i_time]  = integrate_emitted_power(state)

    i_save = something(findfirst(isequal(solution.times[i_time]), solution.output.output_times), 0)
    i_save == 0 && return

    if solution.output.save_surface_temperature
        solution.surface_temperature[:, i_save] .= surface_temperature(state)
    end

    if !isempty(solution.output.subsurface_face_ids)
        for (i, T_sub) in solution.subsurface_temperature
            T_sub[:, i_save] .= state.temperature[:, i]
        end
    end

    if solution.output.save_face_forces
        solution.face_forces[:, i_save] .= state.face_forces
    end

    if !isempty(solution.output.roughness_face_ids)
        for (i, T_rough) in solution.roughness_surface_temperature
            rs = state.roughness_states[state.face_roughness_indices[i]]
            T_rough[:, i_save] .= surface_temperature(rs)
        end
    end
end

"""
    record_timestep!(solution, state, i_time, R)

Record simulation data for a single asteroid at timestep `i_time`, including force and torque
rotated to the inertial frame.

Extends the 3-argument form by additionally recording net thermal force and torque
(when `save_forces`/`save_torques` are `true`) after rotating from the body-fixed frame
to the inertial frame via `R`.

# Arguments
- `solution` : Solution container (`SingleAsteroidThermoPhysicalSolution`)
- `state`    : Current simulation state (`SingleAsteroidThermoPhysicalState`)
- `i_time`   : Index into `solution.times` for the current timestep
- `R`        : Rotation matrix from the body-fixed frame to the inertial frame
"""
function record_timestep!(
    solution ::SingleAsteroidThermoPhysicalSolution,
    state    ::SingleAsteroidThermoPhysicalState,
    i_time   ::Integer,
    R        ::SMatrix{3,3,Float64,9},
)
    record_timestep!(solution, state, i_time)

    i_save = something(findfirst(isequal(solution.times[i_time]), solution.output.output_times), 0)
    i_save == 0 && return

    solution.output.save_forces  && (solution.forces[i_save]  = R * state.force)
    solution.output.save_torques && (solution.torques[i_save] = R * state.torque)
end

"""
    record_timestep!(solution, state, i_time)

Record simulation data for both bodies of a binary system at timestep `i_time`.
Delegates to the single-body form for each body independently.

# Arguments
- `solution` : Solution container (`BinaryAsteroidThermoPhysicalSolution`)
- `state`    : Current simulation state (`BinaryAsteroidThermoPhysicalState`)
- `i_time`   : Index into `solution.primary.times` for the current timestep
"""
function record_timestep!(
    solution ::BinaryAsteroidThermoPhysicalSolution,
    state    ::BinaryAsteroidThermoPhysicalState,
    i_time   ::Integer,
)
    record_timestep!(solution.primary,   state.primary,   i_time)
    record_timestep!(solution.secondary, state.secondary, i_time)
end

"""
    record_timestep!(solution, state, i_time, R₁ᵢ, R₂ᵢ)

Record simulation data for both bodies of a binary system at timestep `i_time`, including
force and torque rotated to the inertial frame for each body.

# Arguments
- `solution` : Solution container (`BinaryAsteroidThermoPhysicalSolution`)
- `state`    : Current simulation state (`BinaryAsteroidThermoPhysicalState`)
- `i_time`   : Index into `solution.primary.times` for the current timestep
- `R₁ᵢ`     : Rotation matrix from the primary body-fixed frame to the inertial frame
- `R₂ᵢ`     : Rotation matrix from the secondary body-fixed frame to the inertial frame
"""
function record_timestep!(
    solution ::BinaryAsteroidThermoPhysicalSolution,
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

    filepath = joinpath(dirpath, "subsurface_temperature.csv")
    CSV.write(filepath, df)
end

# Long format — one row per (time, global face, sub-face) — because roughness models may differ
# in size from face to face, which a wide table with one column per face could not hold.
function _export_roughness_surface_temperature(dirpath, solution::SingleAsteroidThermoPhysicalSolution)
    output_times = solution.output.output_times
    time_col, face_col, sub_col, T_col = Float64[], Int[], Int[], Float64[]
    for i in sort(collect(keys(solution.roughness_surface_temperature)))
        T_rough = solution.roughness_surface_temperature[i]
        n_sub = size(T_rough, 1)
        for (i_save, t) in enumerate(output_times), j in 1:n_sub
            push!(time_col, t); push!(face_col, i); push!(sub_col, j); push!(T_col, T_rough[j, i_save])
        end
    end
    df = DataFrame(time = time_col, face_id = face_col, sub_face_id = sub_col, temperature = T_col)
    filepath = joinpath(dirpath, "roughness_surface_temperature.csv")
    CSV.write(filepath, df)
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

function _export_thermal_net_forces(dirpath, solution::SingleAsteroidThermoPhysicalSolution)
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

"""
    export_solution(dirpath, solution::SingleAsteroidThermoPhysicalSolution)

Export simulation results to CSV files in `dirpath`.
`dirpath` is created automatically if it does not exist.

Files written depend on the `output` specification:
- `diagnostics.csv`            : always (`absorbed_power`, `emitted_power` at all timesteps)
- `surface_temperature.csv`    : when `output.save_surface_temperature = true`
- `subsurface_temperature.csv` : when `output.subsurface_face_ids` is non-empty
- `thermal_face_forces.csv`    : when `output.save_face_forces = true`
- `thermal_net_forces.csv`     : when `output.save_forces = true` or `output.save_torques = true`
- `roughness_surface_temperature.csv` : when `output.roughness_face_ids` is non-empty;
  long format with columns `time`, `face_id`, `sub_face_id`, `temperature`
"""
function export_solution(dirpath, solution::SingleAsteroidThermoPhysicalSolution)
    mkpath(dirpath)
    _export_diagnostics(dirpath, solution)
    solution.output.save_surface_temperature    && _export_surface_temperature(dirpath, solution)
    !isempty(solution.output.subsurface_face_ids) && _export_subsurface_temperature(dirpath, solution)
    solution.output.save_face_forces            && _export_thermal_face_forces(dirpath, solution)
    (solution.output.save_forces || solution.output.save_torques) && _export_thermal_net_forces(dirpath, solution)
    !isempty(solution.output.roughness_face_ids) && _export_roughness_surface_temperature(dirpath, solution)
end


"""
    export_solution(dirpath, solution::BinaryAsteroidThermoPhysicalSolution)

Export results for both bodies to `dirpath/primary/` and `dirpath/secondary/`.
All directories are created automatically if they do not exist.
"""
function export_solution(dirpath, solution::BinaryAsteroidThermoPhysicalSolution)
    dirpath1 = joinpath(dirpath, "primary")
    dirpath2 = joinpath(dirpath, "secondary")
    mkpath(dirpath1)
    mkpath(dirpath2)
    export_solution(dirpath1, solution.primary)
    export_solution(dirpath2, solution.secondary)
end
