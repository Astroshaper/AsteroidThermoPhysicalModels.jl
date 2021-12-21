


"""
"""
function run_YORP(shape, orbit, spin, params_thermo)
    @unpack P, Δt, t_bgn, t_end, Nt, Nz = params_thermo
    
    init_temperature!(shape, orbit, spin, params_thermo)
    
    τ̄ = zeros(3)  # Net YORP torque

    for t in (t_bgn:Δt:t_end)*P
        spin_phase = spin.ω * t
        F☉, r̂☉ = getSolarCondition(orbit, spin, t)
        
        update_flux_sun!(shape, F☉, r̂☉)
        update_flux_scat_single!(shape, params_thermo)
        update_flux_rad_single!(shape, params_thermo)
        
        update_force!(shape, params_thermo)
        # update_force_Rubincam!(shape, params_thermo)
        
        τ̄ .+= body_to_orbit(SVector{3}(shape.torque), spin.γ, spin.ε, spin_phase)
        
        update_temperature!(shape, params_thermo)
        
        @unpack Tz, flux, df = shape.smeshes[1]
        println(Tz[begin], ", ", flux.sun, ", ", flux.scat, ", ", flux.rad, ", ", df)
    end
    τ̄ /= Nt
end


"""
"""
function run_Yarkovsky(shape, orbit, spin, params_thermo)
    @unpack P, Δt, t_bgn, t_end, Nt, Nz = params_thermo
    
    init_temperature!(shape, orbit, spin, params_thermo)
    
    ts = (t_bgn:Δt:t_end)*P
    # Fs = typeof(shape.force)[]
    # Fs = Vector{typeof(shape.force)}(undef, length(ts))
    Fs = [zeros(3) for _ in 1:length(ts)]

    for (i, t) in enumerate(ts)
        spin_phase = spin.ω * t
        F☉, r̂☉ = getSolarCondition(orbit, spin, t)
        
        update_flux_sun!(shape, F☉, r̂☉)
        update_flux_scat_single!(shape, params_thermo)
        update_flux_rad_single!(shape, params_thermo)
        
        update_force!(shape, params_thermo)
        # update_force_Rubincam!(shape, params_thermo)
        
        F = body_to_orbit(SVector{3}(shape.force), spin.γ, spin.ε, spin_phase)
        Fs[i] .= F
        # push!(Fs, F)
        
        update_temperature!(shape, params_thermo)
    end
    Fs
    # τ̄ /= Nt
end

"""
    init_temperature!(shape, orbit, spin, params_thermo)

Initialize temperature distribution of every facet
"""
function init_temperature!(shape, orbit, spin, params_thermo)
    @unpack P, Δt, t_bgn, t_end, Nt, Nz = params_thermo
    
    for facet in shape.facets
        length(facet.Tz) == 0 ? append!(facet.Tz, zeros(Nz)) : facet.Tz .= 0
    end
    length(shape.Tz⁺) == 0 ? append!(shape.Tz⁺, zeros(Nz)) : shape.Tz⁺ .= 0
end


# ****************************************************************
#        Energy flux of sunlight, scattering, and radiation
# ****************************************************************

"""
    update_flux_sun!(shape, F☉, r̂☉)

Update illumination
"""
function update_flux_sun!(shape, F☉, r̂☉)
    for facet in shape.facets
        if isIlluminated(facet, r̂☉, shape)
            facet.flux.sun = F☉ * (facet.normal ⋅ r̂☉)
        else
            facet.flux.sun = 0.
        end
    end
end

"""
    update_flux_scat_single!(shape, params_thermo)

Single scattering of sunlight is considered.
"""
function update_flux_scat_single!(shape, params_thermo)
    @unpack A_B = params_thermo
    
    for facet in shape.facets
        facet.flux.scat = 0.
        for (id, f) in zip(facet.visiblefacets.id, facet.visiblefacets.f)
            facet.flux.scat += f * shape.facets[id].flux.sun
        end
        facet.flux.scat *= A_B
    end
end

"""
    update_flux_scat_mult!(shape, params_thermo)

Multiple scattering of sunlight is considered.
"""
function update_flux_scat_mult!(shape, params_thermo)
    @unpack A_B = params_thermo
    
    # for m in shape.smeshes
    #     m.flux.scat = 0.
    #     for (id, f) in zip(m.visiblefaces.id, m.visiblefaces.f)
    #         m.flux.scat += f * shape.smeshes[id].flux.sun
    #     end
    #     m.flux.scat *= A_B
    # end
end

"""
    update_flux_rad_single!(shape, params_thermo)

Single radiation-reabsorption is considered,
assuming albedo is close to zero at thermal infrared wavelength.
"""
function update_flux_rad_single!(shape, params_thermo)
    @unpack ϵ, A_TH = params_thermo
        
    for facet in shape.facets
        facet.flux.rad = 0.
        for (id, f) in zip(facet.visiblefacets.id, facet.visiblefacets.f)
            T = shape.facets[id].Tz[begin]
            facet.flux.rad += f * T^4
        end
        facet.flux.rad *= ϵ * σ_SB * (1 - A_TH)
    end
end


# ****************************************************************
#                  Photon recoil force / torque
# ****************************************************************

"""
    update_force!(shape, params_thermo)

Update photon recoil force on every facet (df)
"""
function update_force!(shape, params_thermo)
    @unpack A_B, ϵ = params_thermo
    
    for facet in shape.facets
        E = A_B * facet.flux.scat + ϵ * σ_SB * facet.Tz[begin]^4

        @. facet.force = facet.normal
        for vf in facet.visiblefacets
            @. facet.force -= 1.5 * vf.f * vf.d̂
        end
        @. facet.force *= - 2*E*facet.area / (3*c₀)
    end

    shape.force  .= 0
    shape.torque .= 0
    for facet in shape.facets
        r  = SVector{3}(facet.center)
        r̂  = normalize(r)
        df = SVector{3}(facet.force)
        
        shape.force  .+= (r̂ ⋅ df) * r̂  # Photon recoil force
        shape.torque .+= r × df        # Photon recoil torque
    end
end

"""
    update_force_Rubincam!(shape, params_thermo)

Update photon recoil force on every facet (df) based on Rubincam (2000) approximation
"""
function update_force_Rubincam!(shape, params_thermo)
    shape.torque .= 0.
    for facet in shape.facets
        facet.force .= - 2 * facet.flux.sun * facet.area / (3*c₀) .* facet.normal
        shape.torque .+= facet.center × SVector{3}(facet.force)
    end
end


# ****************************************************************
#
# ****************************************************************


"""
    getSolarIrradiation(distance) -> solar_irrad

Calculate the solar irradiation on the body

# Parameter
- `rₕ` : heliocentric distance of the body [m]

# Return
- `F☉` : solar irradiation [W/m^2]
"""
getSolarIrradiation(rₕ) = SOLAR_CONST / (rₕ / AU)^2


"""
    getSolarCondition(orbit, spin, time) -> Φ, r_sun

Get the solar irradition and the direction of the Sun

# Parameters
- `orbit` :
- `spin`  :
- `time`  : [sec]

# Returns
- `F☉` : solar irradiation [W/m^2]
- `r̂☉` : solar direction in the body-fixed frame
"""
function getSolarCondition(orbit, spin, time)
    u = solveKeplerEquation2(orbit, time)
    r = get_r(orbit, u)
    F☉ = getSolarIrradiation(norm(r))

    r̂☉ = normalize(r) * -1  # Shift the origin from the sun to the body

    spin_phase = spin.ω * time
    r̂☉ = orbit_to_body(r̂☉, spin.γ, spin.ε, spin_phase)

    F☉, r̂☉
end


"""

Integrate YORP torque over the entire surface of the shape [Rubincam2000]

# Parameters
- `shape` :
- `r_sun` : direction of the sun in the body-fixed frame (normalized)
- `Φ`     : solar irradiation [W/m^2]

Return
- `τ`

In the future, we need to improve:
- thermal pressure model: Rubincam (2000) model currently implemented.
- local shadow detection
- Tangential YORP effect
"""
function sumTorqueOverSurface(shape, F☉, r̂☉)
    τ = MVector(0., 0., 0.)  # YORP torque

    for mesh in shape.smeshes
        Ψ = mesh.normal ⋅ r̂☉  # cosine of the Sun illumination angle
        if Ψ > 0  # daytime hemisphere of the body
            df = Ψ * mesh.area * mesh.normal  # force on each facet
            dτ = mesh.center × df             # torque on each facet
            τ .+= dτ
        end
    end
    τ *= - 2/3 * F☉ / c₀
    return SVector(τ)
end


function sumTorqueOverSurface_shadowing(shape, F☉, r̂☉)
    τ = MVector(0., 0., 0.)  # YORP torque

    for mesh in shape.smeshes
        Ψ = mesh.normal ⋅ r̂☉  # cosine of the Sun illumination angle
        if Ψ > 0  # daytime hemisphere of the body
            if isIlluminated(mesh, r̂☉, shape.smeshes)
                df = Ψ * mesh.area * mesh.normal  # force on each facet
                dτ = mesh.center × df             # torque on each facet
                τ .+= dτ
            end
        end
    end
    τ *= - 2/3 * F☉ / c₀
    return SVector(τ)
end


"""
    getNetTorque(shape, orbit, spin, times) -> τ̄

Average YORP torque over given time steps
"""
function getNetTorque(shape, orbit, spin, times)

    τ̄ = MVector(0., 0., 0.)  # net YORP torque

    for time in times
        spin_phase = spin.ω * time
        F☉, r̂☉ = getSolarCondition(orbit, spin, time)
        τ = sumTorqueOverSurface(shape, F☉, r̂☉)
        τ = body_to_orbit(τ, spin.γ, spin.ε, spin_phase)

        τ̄ .+= τ
    end
    τ̄ /= length(times)
end

function getNetTorque_shadowing(shape, orbit, spin, times)

    τ̄ = MVector(0., 0., 0.)  # net YORP torque

    for time in times
        spin_phase = spin.ω * time
        F☉, r̂☉ = getSolarCondition(orbit, spin, time)
        τ = sumTorqueOverSurface_shadowing(shape, F☉, r̂☉)
        τ = body_to_orbit(τ, spin.γ, spin.ε, spin_phase)

        τ̄ .+= τ
    end
    τ̄ /= length(times)
end


# ****************************************************************
#
# ****************************************************************


"""
Integrate torque over a whole surface: Pseudo-convex model
"""
function getTotalTorque_convex(shape, F☉, r̂☉)
    τ = MVector(0., 0., 0.)
    for mesh in shape.smeshes
        Ψ = mesh.normal ⋅ r̂☉  # cosine of the Sun illumination angle
        Ψ > 0 && (τ .+= gettorque(mesh, F☉ * Ψ))
    end
    return SVector(τ)
end


"""
Integrate torque over a whole surface: Self-shadowing model
"""
function getTotalTorque_shadow(shape, F☉, r̂☉)
    τ = MVector(0., 0., 0.)
    for mesh in shape.smeshes
        Ψ = mesh.normal ⋅ r̂☉  # cosine of the Sun illumination angle
        Ψ > 0 && isIlluminated(mesh, r̂☉, shape) && (τ .+= gettorque(mesh, F☉ * Ψ))
    end
    return SVector(τ)
end

#####################

"""
Net torque during a revolutional cycle
"""
function getAnnualNetTorque(shape, orbit, spin, times; shadow)
    if shadow
        getTotalTorque = getTotalTorque_shadow
    else
        getTotalTorque = getTotalTorque_convex
    end
    
    getAnnualNetTorque(shape, orbit, spin, times, getTotalTorque)
end


function getAnnualNetTorque(shape, orbit, spin, times, getTotalTorque)
    τ̄ = MVector(0., 0., 0.)

    for time in times
        spin_phase = spin.ω * time
        F☉, r̂☉ = getSolarCondition(orbit, spin, time)
        τ = getTotalTorque(shape, F☉, r̂☉)
        τ = body_to_orbit(τ, spin.γ, spin.ε, spin_phase)

        τ̄ .+= τ
    end
    τ̄ /= length(times)
end






# ****************************************************************
#                            Analysis
# ****************************************************************


getτω(τ, spin) = τ ⋅ spin.ŝ
getτε(τ, spin) = τ ⋅ getSpinUnit1(spin)
getτψ(τ, spin) = τ ⋅ getSpinUnit2(spin)


"""
    torque2rate(τ, spin, C) -> ω̇, ωε̇, ωψ̇

# Parameters
- `τ`    : torque in a body's orbital plane frame
- `spin` :
- `C`    : moment of inertia
"""
function torque2rate(τ, spin, C)
    τ_ω = getτω(τ, spin)
    τ_ε = getτε(τ, spin)
    τ_ψ = getτψ(τ, spin)

    ## [rad/sec/sec]
    ω̇  = τ_ω / C
    ωε̇ = τ_ε / C
    ωψ̇ = τ_ψ / C

    ## [deg/day/day]
    ω̇  = rad2deg(ω̇)  * (3600*24)^2
    ωε̇ = rad2deg(ωε̇) * (3600*24)^2
    ωψ̇ = rad2deg(ωψ̇) * (3600*24)^2

    return ω̇, ωε̇, ωψ̇
end


function torque2rate(τ, spin)
    τ_ω = getτω(τ, spin)
    τ_ε = getτε(τ, spin)
    τ_ψ = getτψ(τ, spin)

    ## [rad/sec/sec]
    ω̇  = τ_ω / shape.I[3, 3]
    ωε̇ = τ_ε / shape.I[3, 3]
    ωψ̇ = τ_ψ / shape.I[3, 3]

    ## [deg/day/day]
    ω̇  = rad2deg(ω̇)  * (3600*24)^2
    ωε̇ = rad2deg(ωε̇) * (3600*24)^2
    ωψ̇ = rad2deg(ωψ̇) * (3600*24)^2

    return ω̇, ωε̇, ωψ̇
end


"""
Caluculate the YORP time scale given a constant acceleration/deceleration

Parameters
- `T_start` [hour]
- `T_end`   [hour]
- `ω̇`       [deg/day/day]
"""
function getTimeScale(T_start, T_end, ω̇)
    T_start *= 3600  # [sec]
    T_end   *= 3600  # [sec]

    ω_start = 2π / T_start  # initial spin rate [rad/sec]
    ω_end   = 2π / T_end    # final spin rate [rad/sec]
    
    Δω = ω_end - ω_start  # [rad/sec]
    Δω = rad2deg(Δω)      # [deg/sec]
    Δω *= 3600*24         # [deg/day]
    
    timescale = Δω / ω̇   # [day]
    timescale /= 365e6   # [Myr]
end


# ****************************************************************
#                   Time variance of torque
# ****************************************************************

"""
    getTorqueVariation(shape, orbit, spin, times) -> τs

Get time variation of YORP torque at each time step

"""
function getTorqueVsTime(shape, orbit, spin, times)
    
    τs = SVector{3,Float64}[]

    for time in times
        spin_phase = spin.ω * time
        F_sun, r̂_sun = getSolarCondition(orbit, spin, time)
        τ = sumTorqueOverSurface(shape, F_sun, r̂_sun)
        τ = body_to_orbit(τ, spin.γ, spin.ε, spin_phase)
        
        push!(τs, τ)
    end
    τs
end

function getTorqueVsTime_shadowing(shape, orbit, spin, times)
    
    τs = SVector{3,Float64}[]

    for time in times
        spin_phase = spin.ω * time
        F_sun, r̂_sun = getSolarCondition(orbit, spin, time)
        τ = sumTorqueOverSurface_shadowing(shape, F_sun, r̂_sun)
        τ = body_to_orbit(τ, spin.γ, spin.ε, spin_phase)
        
        push!(τs, τ)
    end
    τs
end


# ****************************************************************
#              Spatial distribution of net torque
# ****************************************************************


"""
"""
function addTorqueDistribution!(dτ̄s, shape, F_sun, r̂_sun, spin_phase)
    # c₀ = 299792458.0  # speed of light [m/s]

    for i in eachindex(shape.smeshes)
        mesh = shape.smeshes[i]
        Ψ = mesh.normal ⋅ r̂_sun  # cosine of the Sun illumination angle
        if Ψ > 0  # daytime hemisphere of the body
            df = - 2/3 * F_sun * Ψ / c₀ * mesh.area * mesh.normal
            dτ = mesh.center × df
            dτ = rotateZ(dτ, -spin_phase)

            dτ̄s[i] .+= dτ
        end
    end
end

function addTorqueDistribution_shadowing!(dτ̄s, shape, F_sun, r̂_sun, spin_phase)
    # c₀ = 299792458.0  # speed of light [m/s]

    for i in eachindex(shape.smeshes)
        mesh = shape.smeshes[i]
        Ψ = mesh.normal ⋅ r̂_sun  # cosine of the Sun illumination angle
        if Ψ > 0  # daytime hemisphere of the body
            if isIlluminated(mesh, r̂_sun, shape.smeshes)
                df = - 2/3 * F_sun * Ψ / c₀ * mesh.area * mesh.normal
                dτ = mesh.center × df
                dτ = rotateZ(dτ, -spin_phase)  # body's equatorial coordinate

                dτ̄s[i] .+= dτ
            end
        end
    end
end
    

"""
"""
function getNetTorqueDistribution(shape, orbit, spin, times)

    ## Spatial distribution of net torque
    dτ̄s = [MVector(0.,0.,0.) for i in eachindex(shape.smeshes)]

    for time in times
        spin_phase = spin.ω * time
        F_sun, r̂_sun = getSolarCondition(orbit, spin, time)
        addTorqueDistribution!(dτ̄s, shape, F_sun, r̂_sun, spin_phase)
    end
    dτ̄s ./= length(times)

    for dτ̄ in dτ̄s
        rotateX!(dτ̄, -spin.ε)  # body's ecliptic coordinate
        rotateZ!(dτ̄, -spin.γ)  # orbital plane frame
    end
    return dτ̄s
end


function getNetTorqueDistribution_shadowing(shape, orbit, spin, times)

    ## Spatial distribution of net torque
    dτ̄s = [MVector(0.,0.,0.) for i in eachindex(shape.smeshes)]

    for time in times
        spin_phase = spin.ω * time
        F_sun, r̂_sun = getSolarCondition(orbit, spin, time)
        addTorqueDistribution_shadowing!(dτ̄s, shape, F_sun, r̂_sun, spin_phase)
    end
    dτ̄s ./= length(times)

    for dτ̄ in dτ̄s
        rotateX!(dτ̄, -spin.ε)  # body's ecliptic coordinate
        rotateZ!(dτ̄, -spin.γ)  # orbital plane frame
    end
    return dτ̄s
end


