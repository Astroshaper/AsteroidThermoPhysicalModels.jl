#=
thermal_radiation.jl

Tests for thermal radiation calculations.
This test validates:
- Blackbody radiation calculations (Planck's law)
- Stefan-Boltzmann law implementation
- View factor calculations for thermal radiation exchange
- Self-heating effects between asteroid surface elements
- flux_rad being the incident irradiance, independent of the thermal-infrared reflectance
=#

@testset "thermal_radiation" begin

    msg = """\n
    ┏━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━┓
    |                Test: thermal_radiation                 |
    ┗━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━┛
    """
    println(msg)

    ## Calculate the intensity of blackbody radiation of at a wavelength of 6e-7 m and a temperature of 5850 K,
    ## and compare it with the value calculated by Planck.jl.
    ## cf. https://github.com/JuliaAstro/Planck.jl/blob/main/src/Planck.jl
    @test AsteroidThermoPhysicalModels.blackbody_radiance(6e-7, 5850) ≈ 2.583616647617974e13

    ## Check the value of Stefan-Boltzmann law at a temperature of 5850 K
    @test AsteroidThermoPhysicalModels.blackbody_radiance(5850) ≈ AsteroidThermoPhysicalModels.σ_SB * 5850^4

    ######## Thermal radiance from the local terrain model ########

    ## --- Load shape model ---
    path_obj = joinpath("shape", "fractal_v2572_f5000.obj")
    shape = load_shape_obj(path_obj; scale=1, with_face_visibility=true)
    n_face = length(shape.faces)  # Number of faces

    ## --- Ephemerides ---
    P = SPICE.convrt(8, "hours", "seconds")  # Rotation period of the asteroid [s]
    
    ncycles = 2  # Number of cycles to perform TPM
    nsteps_in_cycle = 360  # Number of time steps in one rotation period
    
    ## TPM simulation time duration (ephemerides time)
    et_begin = 0.0  # Start time of TPM
    et_end   = et_begin + P * ncycles  # End time of TPM
    et_range = range(et_begin, et_end; length=nsteps_in_cycle*ncycles+1)

    ## Rotation matrix between the global (asteroid-fixed) frame and the local frame at `(lat, lon)`.
    local_to_global(lat, lon) = RotZYZ(lon, π/2 - lat, π/2)
    global_to_local(lat, lon) = inv(local_to_global(lat, lon))
    
    ## Solar position in the local frame of the roughness map
    lat = deg2rad(0)  # Latitude of the local frame [rad]
    lon = deg2rad(0)  # Longitude of the local frame [rad]

    r☉₀ = [SPICE.convrt(-1, "au", "m"), 0, 0]  # Sun's position at the initial time step in the asteoid-fixed frame

    times = collect(et_range)
    r_sun = [global_to_local(lat, lon) * inv(RotZ(2π * et / P)) * r☉₀ for et in et_range]  # Sun's position in local frame [m]

    ephem = SingleAsteroidEphemerides(times, r_sun)
    
    ## --- Thermal properties ---
    k  = 0.1     # Thermal conductivity [W/m/K]
    ρ  = 1270.0  # Density [kg/m³]
    Cₚ = 600.0   # Heat capacity [J/kg/K]
    
    l = thermal_skin_depth(P, k, ρ, Cₚ)
    Γ = thermal_inertia(k, ρ, Cₚ)

    thermo_params = ThermoParams(
        conductivity    = k,     # Thermal conductivity [W/m/K]
        density         = ρ,     # Density [kg/m³]
        heat_capacity   = Cₚ,    # Heat capacity [J/kg/K]
        reflectance_vis = 0.04,  # Reflectance in visible light [-]
        reflectance_ir  = 0.0,   # Reflectance in thermal infrared [-]
        emissivity      = 1.0,   # Emissivity [-]
    )
    grid_params = GridParams(; z_max=0.6, n_depth=41)

    ## --- Setting of TPM ---
    problem = SingleAsteroidThermoPhysicalProblem(shape, thermo_params, grid_params;
        with_self_shadowing      = true,
        with_self_heating        = true,
        upper_boundary_condition = RadiationBoundaryCondition(),
        lower_boundary_condition = InsulationBoundaryCondition(),
    )

    output_times        = ephem.times[end-nsteps_in_cycle:end]  # Save temperature during the final rotation
    subsurface_face_ids = [49, 340, 648]  # Face indices to save subsurface temperature
    output = SingleAsteroidOutputSpec(output_times; subsurface_face_ids)

    solution = solve(problem, ExplicitEuler();
        ephem               = ephem,
        output              = output,
        initial_temperature = 200.0,
    )

    ## --- Check the thermal radiation from the local terrain model ---
    obs_above = SVector{3, Float64}(0, 0, 1000)  # Observer is just above the local terrain model
    obs_east  = RotY(+π/6) * obs_above           # Observed from 30° east
    obs_west  = RotY(-π/6) * obs_above           # Observed from 30° west

    emissivities = fill(1.0, length(shape.faces))
    temperatures = solution.surface_temperature[:, 181]

    ## Expected values of the thermal radiation
    ## - Observation from 30° east   : 19.89394442347112  [W/m²]
    ## - Observation from just above : 18.01010231251351  [W/m²]
    ## - Observation from 30° west   : 12.406927167050457 [W/m²]
    @test AsteroidThermoPhysicalModels.thermal_radiance(shape, emissivities, temperatures, obs_east)  ≈ 19.89394442347112
    @test AsteroidThermoPhysicalModels.thermal_radiance(shape, emissivities, temperatures, obs_above) ≈ 18.01010231251351
    @test AsteroidThermoPhysicalModels.thermal_radiance(shape, emissivities, temperatures, obs_west)  ≈ 12.406927167050457

    ######## flux_rad is the incident irradiance ########

    @testset "flux_rad is independent of R_ir" begin
        # The thermal radiation a face receives from its neighbours is their emission ε σ T⁴
        # weighted by the view factor. The receiving face's reflectance belongs to the surface
        # boundary condition, not to flux_rad — so flux_rad must not change with R_ir at all.
        # A concave crater makes the exchange non-zero; a temperature gradient makes a wrong
        # face index visible.
        grid_params = GridParams(; z_max=0.1, n_depth=10)
        make_state(R_ir) = begin
            thermo_params = ThermoParams(conductivity=0.1, density=1000.0, heat_capacity=700.0,
                                         reflectance_vis=0.1, reflectance_ir=R_ir, emissivity=0.9)
            shape   = create_shape_crater(0.4, 0.1; Nx=8, Ny=8)
            problem = SingleAsteroidThermoPhysicalProblem(shape, thermo_params, grid_params;
                with_self_shadowing=false, with_self_heating=true)
            state   = AsteroidThermoPhysicalModels._build_single_state(problem, CrankNicolson())
            n_faces = length(shape.faces)
            T₀ = repeat(reshape(range(200.0, 300.0; length=n_faces) |> collect, 1, n_faces), grid_params.n_depth)
            AsteroidThermoPhysicalModels.init_temperature!(state, T₀)
            AsteroidThermoPhysicalModels.update_flux_rad_single!(state)
            state
        end
        state0 = make_state(0.0)
        state3 = make_state(0.3)

        @test any(>(0), state0.flux_rad)
        @test state3.flux_rad == state0.flux_rad

        # Explicit formula: flux_rad[i] = Σⱼ εⱼ σ Tⱼ⁴ Fᵢⱼ over the faces visible from i
        shape = state3.problem.shape
        σ = AsteroidThermoPhysicalModels.σ_SB
        for i in eachindex(shape.faces)
            js = get_visible_face_indices(shape.face_visibility_graph, i)
            fs = get_view_factors(shape.face_visibility_graph, i)
            expected = sum(0.9 * σ * state3.temperature[begin, j]^4 * f for (j, f) in zip(js, fs); init=0.0)
            @test state3.flux_rad[i] ≈ expected
        end

        # The receiving face's reflectance is applied exactly once, where the flux is absorbed
        i = findfirst(>(0), state3.flux_rad)
        @test AsteroidThermoPhysicalModels.absorbed_energy_flux(0.1, 0.3, 0.0, 0.0, state3.flux_rad[i]) ≈ 0.7 * state3.flux_rad[i]
    end

    @testset "mutual heating flux_rad is independent of R_ir" begin
        # Two flat faces facing each other: the secondary sits above the primary, flipped so
        # that its normal points down. Mutual heating must not depend on either R_ir.
        grid_params = GridParams(; z_max=0.1, n_depth=10)
        path_obj = joinpath(@__DIR__, "shape", "single_face.obj")
        make_binary(R_ir) = begin
            thermo_params = ThermoParams(conductivity=0.1, density=1000.0, heat_capacity=700.0,
                                         reflectance_vis=0.1, reflectance_ir=R_ir, emissivity=0.9)
            shape1 = load_shape_obj(path_obj)
            shape2 = load_shape_obj(path_obj)
            problem = BinaryAsteroidThermoPhysicalProblem((shape1, shape2), thermo_params, grid_params;
                with_self_shadowing=false, with_self_heating=false,
                with_mutual_shadowing=false, with_mutual_heating=true)
            state = AsteroidThermoPhysicalModels._build_binary_state(problem, CrankNicolson())
            AsteroidThermoPhysicalModels.init_temperature!(state, 250.0, 300.0)
            r₁₂ = SVector(0.0, 0.0, 1.0)                       # secondary 1 m above the primary
            R₂₁ = SMatrix{3,3}(1.0, 0.0, 0.0, 0.0, -1.0, 0.0, 0.0, 0.0, -1.0)   # flipped
            AsteroidThermoPhysicalModels.mutual_heating!(state, r₁₂, R₂₁)
            state
        end
        state0 = make_binary(0.0)
        state3 = make_binary(0.3)

        @test all(>(0), state0.primary.flux_rad)
        @test all(>(0), state0.secondary.flux_rad)
        @test state3.primary.flux_rad   == state0.primary.flux_rad
        @test state3.secondary.flux_rad == state0.secondary.flux_rad
    end
end
