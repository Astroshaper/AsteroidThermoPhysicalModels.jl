
"""
    blackbody_radiation(λ, T) -> radiance

Accroding to Planck's law, calculate the spectral intensity of blackbody radiation
at wavelength `λ` and temperature `T`.

# Arguments
- `λ` : Wavelength [m]
- `T` : Temperature [K]

# Return
- `radiance` : Spectral radiance [W/m²/m/steradian]

cf. https://github.com/JuliaAstro/Planck.jl/blob/main/src/Planck.jl
"""
blackbody_radiation(λ, T) = 2 * h * c₀^2 / λ^5 / expm1(h * c₀ / (λ * k_B * T))

"""
    blackbody_radiation(T) -> radiance

According to Stefan-Boltzmann law, calculate the total radiance of blackbody radiation
at temperature `T`, integrated over all wavelength.

# Arguments
- `T` : Temperature [K]

# Return
- `radiance` : Total radiance [W/m²]
"""
blackbody_radiation(T) = σ_SB * T^4


"""
    thermal_radiation(shape, emissivities, temperatures, obs) -> radiance

Calculate the radiance from an asteroid based on the shape model and temperature distribution.

# Arguments
- `shape`        : Shape model of an asteroid
- `emissivities` : Emissivity of each facet of the shape model [-]
- `temperatures` : Temperature of each facet of the shape model [K]
- `obs`          : Position vector of the observer in the same coordinate system as `shape` [m]
"""
function thermal_radiation(shape, emissivities, temperatures, obs)
    if length(temperatures) != length(shape.faces)
        throw(ArgumentError("Length of `temperatures` must be equal to the number of faces of the shape model."))
    end

    radiance = 0.0
    for i in eachindex(shape.faces)
        c = shape.face_centers[i]
        n̂ = shape.face_normals[i]
        a = shape.face_areas[i]

        ε = emissivities[i]
        T = temperatures[i]
        
        ## Direction and distance from the facet to the observer
        d̂ = normalize(obs - c)
        d = norm(obs - c)

        cosθ = n̂ ⋅ d̂
        cosθ < 0 && continue  # Ignore the facet if the observer is on the back side

        ## TO DO: Ray-trace for each facet

        radiance += ε * σ_SB * T^4 * a * cosθ / (π * d^2)  # π, 2π, or 4π?
    end
    return radiance
end

