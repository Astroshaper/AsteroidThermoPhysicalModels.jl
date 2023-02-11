# Test of calculating the Keplerian motion of asteroid Ryugu.
# These functions are used in https://github.com/MasanoriKanamaru/Astroshaper-examples/blob/main/TPM_Kanamaru2021/TPM_Kanamaru2021.ipynb

@testset "Keplerian_motion" begin
    orbit = Astroshaper.OrbitalElements(Astroshaper.RYUGU)
    spin = Astroshaper.SpinParams(Astroshaper.RYUGU, orbit)

    Δt = spin.P / 72
    times = collect(0:Δt:orbit.P)

    for time in times
        spin_phase = spin.ω * time

        u = Astroshaper.solveKeplerEquation2(orbit, time)
        r = Astroshaper.get_r(orbit, u)
        F☉ = Astroshaper.SOLAR_CONST / (norm(r) / Astroshaper.AU)^2
        
        r̂☉ = normalize(r) * -1  # Shift the origin from the sun to the body
        r̂☉ = Astroshaper.orbit_to_body(r̂☉, spin.γ, spin.ε, spin_phase)  # Sun's position in the asteroid-fixed frame
    end
end