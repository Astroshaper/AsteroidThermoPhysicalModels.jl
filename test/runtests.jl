using Astroshaper
using Test

@testset "Astroshaper.jl" begin
    shapepath = "./ryugu_test.obj"
    shape = setShapeModel(shapepath; scale=1000)
    
    params_orbit = Dict()
    params_orbit[:a] = 1.18956373  # semi-mojor axis [AU]
    params_orbit[:e] = 0.19027921  # eccentricity
    params_orbit[:I] = 5.8840222   # inclination [deg]
    params_orbit[:Ω] = 251.589203  # longitude of the ascending node [deg]
    params_orbit[:ω] = 211.435963  # argument of periapsis [deg]
    params_orbit[:Φ] = 21.9353799  # mean anomaly [deg]
    params_orbit[:μ] = GM☉ + 30.0
    
    orbit = OrbitalElements(params_orbit)
    
    params_spin = Dict()
    params_spin[:α] = 96.4
    params_spin[:δ] = -66.4
    params_spin[:T] = 7.63262
    
    spin = setSpinParams(params_spin, orbit)
    
    dt = spin.T / 72
    times = Vector(0:dt:orbit.T)
    
    @time τ̄ = getNetTorque(shape, orbit, spin, times)
    
    @test true
end

