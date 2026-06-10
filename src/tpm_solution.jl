#=
tpm_solution.jl

Solution types and export functions for thermophysical simulations.
=#

# ╔═══════════════════════════════════════════════════════════════════╗
# ║                     Output data format                            ║
# ╚═══════════════════════════════════════════════════════════════════╝

"""
    struct SingleAsteroidThermoPhysicalSolution

Output data format for `SingleAsteroidThermoPhysicalState`

# Fields
## Saved at all time steps
- `times`  : Timesteps, given the same vector as `ephem.time` [s]
- `E_in`   : Input energy per second on the whole surface [W]
- `E_out`  : Output energy per second from the whole surface [W]
- `force`  : Thermal force on the asteroid [N]
- `torque` : Thermal torque on the asteroid [N ⋅ m]

## Saved only at the time steps desired by the user
- `times_to_save`          : Timesteps to save temperature and thermal force on every face [s]
- `depth_nodes`            : Depths of the calculation nodes for 1-D heat conduction [m], a vector of size `n_depth`
- `surface_temperature`    : Surface temperature [K], a matrix in size of `(n_face, n_time)`.
    - `n_face` : Number of faces
    - `n_time` : Number of time steps to save surface temperature
- `subsurface_temperature` : Temperature [K] as a function of depth [m] and time [s], `Dict` with face ID as key and a matrix `(n_depth, n_time)` as an entry.
    - `n_depth` : The number of the depth nodes
    - `n_time` : The number of time steps to save temperature
- `face_forces`            : Thermal force on every face of the shape model [N], a matrix in size of `(n_face, n_time)`.
    - `n_face` : Number of faces
    - `n_time` : Number of time steps to save surface temperature
"""
struct SingleAsteroidThermoPhysicalSolution
    times  ::Vector{Float64}
    E_in   ::Vector{Float64}
    E_out  ::Vector{Float64}
    force  ::Vector{SVector{3, Float64}}
    torque ::Vector{SVector{3, Float64}}

    times_to_save          ::Vector{Float64}
    depth_nodes            ::Vector{Float64}
    surface_temperature    ::Matrix{Float64}
    subsurface_temperature ::Dict{Int, Matrix{Float64}}
    face_forces            ::Matrix{SVector{3, Float64}}
end


"""
Outer constructor of `SingleAsteroidThermoPhysicalSolution`

# Arguments
- `state`         : Thermophysical simulation state for a single asteroid
- `ephem`         : Ephemerides
- `times_to_save` : Timesteps to save temperature
- `face_ID`       : Face indices to save subsurface temperature
"""
function SingleAsteroidThermoPhysicalSolution(state::SingleAsteroidThermoPhysicalState, ephem, times_to_save::Vector{Float64}, face_ID::Vector{Int})
    n_step = length(ephem.time)             # Number of time steps
    n_step_to_save = length(times_to_save)  # Number of time steps to save temperature
    n_face = length(state.problem.shape.faces)       # Number of faces of the shape model

    E_in   = zeros(n_step)
    E_out  = zeros(n_step)
    force  = zeros(SVector{3, Float64}, n_step)
    torque = zeros(SVector{3, Float64}, n_step)

    depth_nodes = state.problem.thermo_params.Δz * (0:state.problem.thermo_params.n_depth-1)
    surface_temperature = zeros(n_face, n_step_to_save)
    subsurface_temperature = Dict{Int,Matrix{Float64}}(
        i => zeros(state.problem.thermo_params.n_depth, n_step_to_save) for i in face_ID
    )
    face_forces = zeros(SVector{3, Float64}, n_face, n_step_to_save)

    return SingleAsteroidThermoPhysicalSolution(
        ephem.time,
        E_in,
        E_out,
        force,
        torque,
        times_to_save,
        depth_nodes,
        surface_temperature,
        subsurface_temperature,
        face_forces,
    )
end


"""
    struct BinaryAsteroidThermoPhysicalSolution

Output data format for `BinaryAsteroidThermoPhysicalState`

# Fields
- `primary`   : Solution for the primary
- `secondary` : Solution for the secondary
"""
struct BinaryAsteroidThermoPhysicalSolution
    primary  ::SingleAsteroidThermoPhysicalSolution
    secondary::SingleAsteroidThermoPhysicalSolution
end


"""
Outer constructor of `BinaryAsteroidThermoPhysicalSolution`

# Arguments
- `state`         : Thermophysical simulation state for a binary asteroid
- `ephem`         : Ephemerides
- `times_to_save` : Timesteps to save temperature (Common to both the primary and the secondary)
- `face_ID_pri`   : Face indices to save subsurface temperature of the primary
- `face_ID_sec`   : Face indices to save subsurface temperature of the secondary
"""
function BinaryAsteroidThermoPhysicalSolution(state::BinaryAsteroidThermoPhysicalState, ephem, times_to_save::Vector{Float64}, face_ID_pri::Vector{Int}, face_ID_sec::Vector{Int})
    result_primary   = SingleAsteroidThermoPhysicalSolution(state.primary,   ephem, times_to_save, face_ID_pri)
    result_secondary = SingleAsteroidThermoPhysicalSolution(state.secondary, ephem, times_to_save, face_ID_sec)

    return BinaryAsteroidThermoPhysicalSolution(result_primary, result_secondary)
end


"""
    record_timestep!(solution::SingleAsteroidThermoPhysicalSolution, state::SingleAsteroidThermoPhysicalState, i_time::Integer)

Save the results of TPM at the time step `i_time` to `solution`.

# Arguments
- `solution` : Solution object for a single asteroid
- `state`    : Simulation state for a single asteroid
- `i_time`   : Time step to save data
"""
function record_timestep!(solution::SingleAsteroidThermoPhysicalSolution, state::SingleAsteroidThermoPhysicalState, i_time::Integer)
    solution.E_in[i_time]   = energy_in(state)
    solution.E_out[i_time]  = energy_out(state)
    solution.force[i_time]  = state.force
    solution.torque[i_time] = state.torque

    t = solution.times[i_time]

    if t in solution.times_to_save
        i_time_save = findfirst(isequal(t), solution.times_to_save)

        solution.surface_temperature[:, i_time_save] .= surface_temperature(state)

        for (i, temperature) in solution.subsurface_temperature
            temperature[:, i_time_save] .= state.temperature[:, i]
        end

        solution.face_forces[:, i_time_save] .= state.face_forces
    end
end


"""
    record_timestep!(solution::BinaryAsteroidThermoPhysicalSolution, state::BinaryAsteroidThermoPhysicalState, i_time::Integer)

Save the results of TPM at the time step `i_time` to `solution`.

# Arguments
- `solution` : Solution object for a binary asteroid
- `state`    : Simulation state for a binary asteroid
- `i_time`   : Time step
"""
function record_timestep!(solution::BinaryAsteroidThermoPhysicalSolution, state::BinaryAsteroidThermoPhysicalState, i_time::Integer)
    record_timestep!(solution.primary,   state.primary,   i_time)
    record_timestep!(solution.secondary, state.secondary, i_time)
end


"""
    export_solution(dirpath, result::SingleAsteroidThermoPhysicalSolution)

Export the result of `SingleAsteroidThermoPhysicalState` to CSV files.
The output files are saved in the following directory structure:

    dirpath
    ├── physical_quantities.csv
    ├── subsurface_temperature.csv
    ├── surface_temperature.csv
    └── thermal_force.csv

# Arguments
- `dirpath` : Path to the directory to save CSV files.
- `result`  : Output data format for `SingleAsteroidThermoPhysicalState`
"""
function export_solution(dirpath, result::SingleAsteroidThermoPhysicalSolution)

    df = DataFrame()
    df.time     = result.times
    df.E_in     = result.E_in
    df.E_out    = result.E_out
    df.force_x  = [f[1] for f in result.force]
    df.force_y  = [f[2] for f in result.force]
    df.force_z  = [f[3] for f in result.force]
    df.torque_x = [τ[1] for τ in result.torque]
    df.torque_y = [τ[2] for τ in result.torque]
    df.torque_z = [τ[3] for τ in result.torque]

    CSV.write(joinpath(dirpath, "physical_quantities.csv"), df)

    ##= Surface temperature =##
    filepath = joinpath(dirpath, "surface_temperature.csv")
    df = hcat(
        DataFrame(time = result.times_to_save),
        DataFrame(
            result.surface_temperature',
            ["face_$(i)" for i in 1:size(result.surface_temperature, 1)]
        ),
    )

    CSV.write(filepath, df)

    ##= Subsurface temperature =##
    filepath = joinpath(dirpath, "subsurface_temperature.csv")

    nrows = length(result.depth_nodes) * length(result.times_to_save)
    df = DataFrame(
        time  = reshape([t for _ in result.depth_nodes, t in result.times_to_save], nrows),
        depth = reshape([d for d in result.depth_nodes, _ in result.times_to_save], nrows),
    )

    # Add a column for each face
    for (i, subsurface_temperature) in collect(result.subsurface_temperature)
        df[:, "face_$(i)"] =
            reshape(subsurface_temperature, length(subsurface_temperature))
    end

    # Sort the columns by the face ID
    keys_sorted = sort(names(df[:, 3:end]), by=x->parse(Int, replace(x, r"[^0-9]" => "")))
    df = df[:, ["time", "depth", keys_sorted...]]

    CSV.write(filepath, df)

    ##= Thermal force on every face of the shape model =##
    filepath = joinpath(dirpath, "thermal_force.csv")

    n_face = size(result.face_forces, 1)  # Number of faces of the shape model
    n_step = size(result.face_forces, 2)  # Number of time steps to save temperature
    nrows = n_face * n_step

    df = DataFrame(
        time = reshape([t for _ in 1:n_face, t in result.times_to_save], nrows),
        face = reshape([i for i in 1:n_face, _ in result.times_to_save], nrows),
    )
    df.x = reshape([f[1] for f in result.face_forces], nrows)  # x-component of the thermal force
    df.y = reshape([f[2] for f in result.face_forces], nrows)  # y-component of the thermal force
    df.z = reshape([f[3] for f in result.face_forces], nrows)  # z-component of the thermal force

    CSV.write(filepath, df)
end


"""
    export_solution(dirpath, result::BinaryAsteroidThermoPhysicalSolution)

Export the result of `BinaryAsteroidThermoPhysicalState` to CSV files.
The output files are saved in the following directory structure:

    dirpath
    ├── primary
    │   ├── physical_quantities.csv
    │   ├── subsurface_temperature.csv
    │   ├── surface_temperature.csv
    │   └── thermal_force.csv
    └── secondary
        ├── physical_quantities.csv
        ├── subsurface_temperature.csv
        ├── surface_temperature.csv
        └── thermal_force.csv

# Arguments
- `dirpath`  : Path to the directory to save CSV files.
- `solution` : Solution object for a binary asteroid
"""
function export_solution(dirpath, solution::BinaryAsteroidThermoPhysicalSolution)
    dirpath1 = joinpath(dirpath, "primary")
    dirpath2 = joinpath(dirpath, "secondary")

    mkpath(dirpath1)
    mkpath(dirpath2)

    export_solution(dirpath1, solution.primary)
    export_solution(dirpath2, solution.secondary)
end
