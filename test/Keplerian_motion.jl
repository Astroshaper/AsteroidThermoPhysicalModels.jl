"""
Test of calculating the Keplerian motion of asteroid Ryugu.

These functions are used in https://github.com/MasanoriKanamaru/Astroshaper-examples/blob/main/TPM_Kanamaru2021/TPM_Kanamaru2021.ipynb
"""
@testset "Keplerian_motion" begin

    orbit = Astroshaper.OrbitalElements(Astroshaper.RYUGU)
    spin = Astroshaper.SpinParams(Astroshaper.RYUGU, orbit)

    Δt = spin.P / 72
    times = collect(0:Δt:orbit.P)
    
end