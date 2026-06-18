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

Encapsulates which timesteps and face indices to record in detail,
replacing the individual `times_to_save` and `face_ID` keyword arguments of `solve`.

# Fields
- `times_to_save` : Timesteps at which detailed temperature data are saved [s]; must be a subset of `ephem.times`
- `face_ID`       : Face indices for which to save subsurface temperature profiles
"""
struct SingleAsteroidOutputSpec
    times_to_save ::Vector{Float64}
    face_ID       ::Vector{Int}
end

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

# Validate that all times_to_save are present in ephem.times.
function _validate_output_spec(output::SingleAsteroidOutputSpec, ephem)
    for t in output.times_to_save
        findfirst(isequal(t), ephem.times) === nothing &&
            throw(ArgumentError("times_to_save contains $t which is not in ephem.times"))
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
- `times`   : Timesteps [s]
- `E_in`    : Input power on the whole surface [W]
- `E_out`   : Output power from the whole surface [W]
- `forces`  : Net thermal force in the inertial frame [N], or `nothing`
- `torques` : Net thermal torque in the inertial frame [N⋅m], or `nothing`

## Saved only at timesteps requested via `times_to_save`
- `times_to_save`          : Timesteps at which detailed data are saved [s]
- `depth_nodes`            : Depth of each calculation node [m], size `(n_depth,)`
- `surface_temperature`    : Surface temperature [K], size `(n_face, n_time)`
- `subsurface_temperature` : Subsurface temperature [K] by face ID; each entry is `(n_depth, n_time)`
- `face_forces`            : Per-face thermal force in the body-fixed frame [N], size `(n_face, n_time)`
"""
struct SingleAsteroidThermoPhysicalSolution{F <: Union{Nothing, AbstractVector}}
    times   ::Vector{Float64}
    E_in    ::Vector{Float64}
    E_out   ::Vector{Float64}
    forces  ::F
    torques ::F

    times_to_save          ::Vector{Float64}
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
    n_step         = length(times)
    n_step_to_save = length(output.times_to_save)
    n_face         = length(state.problem.shape.faces)
    n_depth        = state.problem.thermo_params.n_depth

    E_in  = zeros(n_step)
    E_out = zeros(n_step)
    depth_nodes            = state.problem.thermo_params.Δz * collect(0:n_depth-1)
    surface_temperature    = zeros(n_face, n_step_to_save)
    subsurface_temperature = Dict{Int,Matrix{Float64}}(
        i => zeros(n_depth, n_step_to_save) for i in output.face_ID
    )
    face_forces = zeros(SVector{3, Float64}, n_face, n_step_to_save)

    SingleAsteroidThermoPhysicalSolution{Nothing}(
        times, E_in, E_out, nothing, nothing,
        output.times_to_save, depth_nodes, surface_temperature, subsurface_temperature, face_forces,
    )
end

# Allocate a SingleAsteroidThermoPhysicalSolution{Vector{SVector{3,Float64}}} (with force/torque).
function _alloc_solution_with_force(
    state  ::SingleAsteroidThermoPhysicalState,
    times  ::Vector{Float64},
    output ::SingleAsteroidOutputSpec,
)
    n_step         = length(times)
    n_step_to_save = length(output.times_to_save)
    n_face         = length(state.problem.shape.faces)
    n_depth        = state.problem.thermo_params.n_depth

    E_in    = zeros(n_step)
    E_out   = zeros(n_step)
    forces  = zeros(SVector{3, Float64}, n_step)
    torques = zeros(SVector{3, Float64}, n_step)
    depth_nodes            = state.problem.thermo_params.Δz * collect(0:n_depth-1)
    surface_temperature    = zeros(n_face, n_step_to_save)
    subsurface_temperature = Dict{Int,Matrix{Float64}}(
        i => zeros(n_depth, n_step_to_save) for i in output.face_ID
    )
    face_forces = zeros(SVector{3, Float64}, n_face, n_step_to_save)

    SingleAsteroidThermoPhysicalSolution{Vector{SVector{3,Float64}}}(
        times, E_in, E_out, forces, torques,
        output.times_to_save, depth_nodes, surface_temperature, subsurface_temperature, face_forces,
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

# Shared helper: record temperature data at times_to_save.
function _record_temperature!(solution::SingleAsteroidThermoPhysicalSolution, state::SingleAsteroidThermoPhysicalState, i_time::Integer)
    t = solution.times[i_time]
    t in solution.times_to_save || return

    i_save = findfirst(isequal(t), solution.times_to_save)
    solution.surface_temperature[:, i_save] .= surface_temperature(state)
    for (i, T_sub) in solution.subsurface_temperature
        T_sub[:, i_save] .= state.temperature[:, i]
    end
    solution.face_forces[:, i_save] .= state.face_forces
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
    solution.E_in[i_time]  = energy_in(state)
    solution.E_out[i_time] = energy_out(state)
    _record_temperature!(solution, state, i_time)
end

"""
    record_timestep!(solution::SingleAsteroidThermoPhysicalSolution{<:AbstractVector}, state, i_time, R)

Record energy balance, temperature data, and force/torque rotated to the inertial frame via `R`.
"""
function record_timestep!(
    solution ::SingleAsteroidThermoPhysicalSolution{<:AbstractVector},
    state    ::SingleAsteroidThermoPhysicalState,
    i_time   ::Integer,
    R        ::SMatrix{3,3,Float64,9},
)
    solution.E_in[i_time]    = energy_in(state)
    solution.E_out[i_time]   = energy_out(state)
    solution.forces[i_time]  = R * state.force
    solution.torques[i_time] = R * state.torque
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

# Shared helpers for non-force CSV outputs.
function _export_surface_temperature(dirpath, result::SingleAsteroidThermoPhysicalSolution)
    df = hcat(
        DataFrame(time = result.times_to_save),
        DataFrame(
            result.surface_temperature',
            ["face_$(i)" for i in 1:size(result.surface_temperature, 1)],
        ),
    )
    CSV.write(joinpath(dirpath, "surface_temperature.csv"), df)
end

function _export_subsurface_temperature(dirpath, result::SingleAsteroidThermoPhysicalSolution)
    nrows = length(result.depth_nodes) * length(result.times_to_save)
    df = DataFrame(
        time  = reshape([t for _ in result.depth_nodes, t in result.times_to_save], nrows),
        depth = reshape([d for d in result.depth_nodes, _ in result.times_to_save], nrows),
    )
    for (i, T_sub) in collect(result.subsurface_temperature)
        df[:, "face_$(i)"] = reshape(T_sub, length(T_sub))
    end
    keys_sorted = sort(names(df[:, 3:end]), by=x->parse(Int, replace(x, r"[^0-9]" => "")))
    df = df[:, ["time", "depth", keys_sorted...]]
    CSV.write(joinpath(dirpath, "subsurface_temperature.csv"), df)
end

function _export_thermal_force(dirpath, result::SingleAsteroidThermoPhysicalSolution)
    n_face = size(result.face_forces, 1)
    n_step = size(result.face_forces, 2)
    nrows  = n_face * n_step
    df = DataFrame(
        time = reshape([t for _ in 1:n_face, t in result.times_to_save], nrows),
        face = reshape([i for i in 1:n_face, _ in result.times_to_save], nrows),
    )
    df.x = reshape([f[1] for f in result.face_forces], nrows)
    df.y = reshape([f[2] for f in result.face_forces], nrows)
    df.z = reshape([f[3] for f in result.face_forces], nrows)
    CSV.write(joinpath(dirpath, "thermal_force.csv"), df)
end


"""
    export_solution(dirpath, result::SingleAsteroidThermoPhysicalSolution{Nothing})

Export temperature results to CSV (no force/torque columns).
"""
function export_solution(dirpath, result::SingleAsteroidThermoPhysicalSolution{Nothing})
    df = DataFrame(
        time  = result.times,
        E_in  = result.E_in,
        E_out = result.E_out,
    )
    CSV.write(joinpath(dirpath, "physical_quantities.csv"), df)

    _export_surface_temperature(dirpath, result)
    _export_subsurface_temperature(dirpath, result)
    _export_thermal_force(dirpath, result)
end


"""
    export_solution(dirpath, result::SingleAsteroidThermoPhysicalSolution{<:AbstractVector})

Export temperature and force/torque results to CSV.
"""
function export_solution(dirpath, result::SingleAsteroidThermoPhysicalSolution{<:AbstractVector})
    df = DataFrame(
        time     = result.times,
        E_in     = result.E_in,
        E_out    = result.E_out,
        force_x  = [f[1] for f in result.forces],
        force_y  = [f[2] for f in result.forces],
        force_z  = [f[3] for f in result.forces],
        torque_x = [τ[1] for τ in result.torques],
        torque_y = [τ[2] for τ in result.torques],
        torque_z = [τ[3] for τ in result.torques],
    )
    CSV.write(joinpath(dirpath, "physical_quantities.csv"), df)

    _export_surface_temperature(dirpath, result)
    _export_subsurface_temperature(dirpath, result)
    _export_thermal_force(dirpath, result)
end


"""
    export_solution(dirpath, solution::BinaryAsteroidThermoPhysicalSolution)

Export results for both bodies to `dirpath/primary/` and `dirpath/secondary/`.
"""
function export_solution(dirpath, solution::BinaryAsteroidThermoPhysicalSolution)
    dirpath_pri = joinpath(dirpath, "primary")
    dirpath_sec = joinpath(dirpath, "secondary")
    mkpath(dirpath_pri)
    mkpath(dirpath_sec)
    export_solution(dirpath_pri, solution.primary)
    export_solution(dirpath_sec, solution.secondary)
end
