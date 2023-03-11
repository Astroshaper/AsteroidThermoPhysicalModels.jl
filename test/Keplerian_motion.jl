# Test of calculating the Keplerian motion of asteroid Ryugu.
# These functions are used in https://github.com/Astroshaper/Astroshaper-examples/blob/main/TPM_Kanamaru2021/TPM_Kanamaru2021.ipynb

@testset "Keplerian_motion" begin
    orbit = AsteroidThermoPhysicalModels.OrbitalElements(AsteroidThermoPhysicalModels.RYUGU)
    spin = AsteroidThermoPhysicalModels.SpinParams(AsteroidThermoPhysicalModels.RYUGU, orbit)

    Δt = spin.P / 72
    times = collect(0:Δt:orbit.P)

    for time in times
        spin_phase = spin.ω * time

        u = AsteroidThermoPhysicalModels.solveKeplerEquation2(orbit, time)
        r = AsteroidThermoPhysicalModels.get_r(orbit, u)
        F☉ = AsteroidThermoPhysicalModels.SOLAR_CONST / (norm(r) / AsteroidThermoPhysicalModels.AU)^2
        
        r̂☉ = normalize(r) * -1  # Shift the origin from the sun to the body
        r̂☉ = AsteroidThermoPhysicalModels.orbit_to_body(r̂☉, spin.γ, spin.ε, spin_phase)  # Sun's position in the asteroid-fixed frame
    end
end