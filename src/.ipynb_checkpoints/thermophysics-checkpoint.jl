

# ****************************************************************
#                      1D heat conduction
# ****************************************************************


"""
`l_2π` : Thermal skin depth
"""
getThermalSkinDepth(P, k, ρ, Cₚ) = √(4π * P * k / (ρ * Cₚ))


"""
`Γ` : Thermal inertia
"""
getThermalInertia(k, ρ, Cₚ) = √(k * ρ * Cₚ)


"""
- `Tⱼ`   : Temperatures
- `Tⱼ₊₁` : Temperatures at the next time step

i : index of depth
j : index of time step
"""
function forward!(Tⱼ, Tⱼ₊₁, Δτ, Δz, F_total, Γ, P, ϵ)
    λ = 1/4π * (Δτ/Δz^2)

    Tⱼ₊₁[begin+1:end-1] .= (1-2λ)*Tⱼ[begin+1:end-1] + λ*(Tⱼ[begin+2:end] + Tⱼ[begin:end-2])
    # for i in 2:length(Tⱼ)-1
    #     Tⱼ₊₁[i] = (1-2λ)*Tⱼ[i] + λ*(Tⱼ[i+1] + Tⱼ[i-1])
    # end

    Tⱼ₊₁[end] = Tⱼ₊₁[end-1]  # Internal boundary condition: Insulation
    updateSurfaceTemperature!(F_total, Γ, P, Δz, ϵ, Tⱼ₊₁)  # Sourface boundary condition

    Tⱼ .= Tⱼ₊₁
end


"""
Solve Newton's method to get surface temperature 
"""
function updateSurfaceTemperature!(F_total, Γ, P, Δz, ϵ, T)
    for _ in 1:20
        T_pri = T[begin]

        T[begin] -= f(F_total, Γ, P, Δz, ϵ, T) / f_deriv(Γ, P, Δz, ϵ, T)
        err = abs(1 - T_pri / T[begin])
        # println(i,  " : ", err)
        err < 1e-10 && return
    end
end

f(F_total, Γ, P, Δz, ϵ, T) = F_total + Γ / √(4π * P) * (T[begin+1] - T[begin]) / Δz - ϵ*σ*T[begin]^4

f_deriv(Γ, P, Δz, ϵ, T) = - Γ / √(4π * P) / Δz - 4*ϵ*σ*T[begin]^3





# ****************************************************************
#
# ****************************************************************


"""
Intensity of radiation at a wavelength λ and tempertature T
according to the Planck function
"""
function getintensity(λ, T)
    h = 6.62607015e-34  # Planck constant [J⋅s]
    k = 1.380649e-23    # Boltzmann's constant [J/K]

    I = 2 * h * c^2 / λ^5 / (exp(h * c₀ / (λ * k * T)) - 1)
end


ν2λ(ν) = c₀ / ν
λ2ν(λ) = c₀ / λ

