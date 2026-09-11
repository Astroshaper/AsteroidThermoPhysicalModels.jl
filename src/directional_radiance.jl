#=
directional_radiance.jl

Direction-dependent thermal radiance and brightness temperature of the facets of a shape
model, computed as a post-processing step from a solution. A facet with a roughness model
radiates anisotropically (thermal-infrared beaming); this is where that anisotropy is
evaluated for a given observer direction.
=#

# Blackbody radiance per unit area and steradian of a Lambertian emitter: total (σT⁴/π) when
# `λ === nothing`, spectral (Planck at wavelength `λ`) otherwise.
_lambert_radiance(T, ::Nothing) = σ_SB * T^4 / π
_lambert_radiance(T, λ::Real)   = blackbody_radiance(λ, T) / π

# Blackbody temperature whose Lambertian radiance is `L`; inverse of `_lambert_radiance`.
_brightness_temperature(L, ::Nothing) = (π * L / σ_SB)^(1/4)
function _brightness_temperature(L, λ::Real)
    # π L = 2 h c² / λ⁵ / expm1(h c / (λ k T))  ⇒  T = h c / (λ k log1p(2 h c² / (λ⁵ π L)))
    L ≤ 0 && return L == 0 ? 0.0 : NaN
    h * c₀ / (λ * k_B * log1p(2 * h * c₀^2 / (λ^5 * π * L)))
end


"""
    roughness_radiance(shape::ShapeModel, i, T_sub, ε, d̂; λ=nothing) -> L

Thermal radiance of facet `i` of `shape`, whose roughness model has the sub-facet
surface temperatures `T_sub`, towards the observer direction `d̂` given in the body-fixed frame.

The roughness model is a representative patch of the facet's surface, so the radiance is the
emission of its visible sub-facets per unit *projected* area of the patch:

```math
L_i(\\hat{\\mathbf d}) = \\frac{1}{A_\\mathrm{proj}\\,(\\hat{\\mathbf z}\\cdot\\hat{\\mathbf d}_\\mathrm{local})}
\\sum_j V_j(\\hat{\\mathbf d}_\\mathrm{local})\\,(\\hat{\\mathbf n}_j\\cdot\\hat{\\mathbf d}_\\mathrm{local})^+\\,a_j\\,\\frac{\\varepsilon\\,B(T_j)}{\\pi}
```

where ``\\hat{\\mathbf d}_\\mathrm{local}`` is `d̂` rotated into the local frame of the facet,
``V_j`` is 1 when sub-facet `j` is visible from that direction (not hidden by the crater walls),
``A_\\mathrm{proj} = \\sum_j a_j (\\hat{\\mathbf n}_j \\cdot \\hat{\\mathbf z})`` is the projected area of
the patch, and ``B(T) = \\sigma T^4`` (total) or the Planck function at wavelength `λ` (spectral).
For an isothermal patch without shadowing this reduces to the Lambertian ``\\varepsilon B(T)/\\pi``;
a sunlit crater whose hot wall faces the observer radiates more than that — thermal-infrared
beaming.

# Arguments
- `shape` : Shape model with surface roughness; facet `i` must carry a roughness model
- `i`     : Global facet index
- `T_sub` : Surface temperature of each sub-facet of the roughness model [K]
- `ε`     : Emissivity of the facet (grey: independent of wavelength)
- `d̂`     : Direction from the facet to the observer in the body-fixed frame (normalised internally)

# Keyword Arguments
- `λ` : Wavelength [m] for the spectral radiance; `nothing` (default) for the total radiance

# Returns
- `L` : Radiance [W/m²/sr], or spectral radiance [W/m²/m/sr] when `λ` is given. `NaN` when the
  facet is seen from behind (``\\hat{\\mathbf z}\\cdot\\hat{\\mathbf d}_\\mathrm{local} \\le 0``)

# Notes
- Sub-facet visibility is evaluated with `update_illumination!` of the roughness model with the
  observer direction in place of the Sun: being lit from a direction and being visible from it
  are the same test. The roughness model's visibility graph and maximum elevations are built on
  demand if missing.
- The patch has no neighbours, so at grazing angles rays that would be blocked by the next
  patch are not; the representative-patch picture assumes patches much smaller than the facet.
"""
function roughness_radiance(
    shape::ShapeModel, i::Integer, T_sub::AbstractVector{<:Real}, ε::Real,
    d̂::StaticVector{3}; λ::Union{Nothing, Real} = nothing,
)
    patch = get_roughness_model(shape, i)::ShapeModel
    length(T_sub) == length(patch.faces) || throw(ArgumentError(
        "T_sub has $(length(T_sub)) entries but the roughness model of facet $i has $(length(patch.faces)) facets"
    ))
    d̂_local = transform_physical_vector_global_to_local(shape, i, normalize(d̂))
    cosθ_view = d̂_local[3]
    cosθ_view ≤ 0 && return NaN

    _prepare_self_shadowing!(patch)
    visible = Vector{Bool}(undef, length(patch.faces))
    update_illumination!(visible, patch, d̂_local; with_self_shadowing=true)

    return _directional_emission(patch, visible, d̂_local, n -> ε * _lambert_radiance(T_sub[n], λ))
end


"""
    directional_radiance(problem, solution, i_save, d̂; λ=nothing) -> L

Thermal radiance of every global facet towards the observer direction `d̂` (body-fixed frame)
at output time `solution.output.output_times[i_save]`.

Facets whose roughness-model surface temperatures were recorded
(`output.roughness_face_ids`) radiate
anisotropically according to [`roughness_radiance`](@ref); every other facet radiates as a
smooth Lambertian surface, ``\\varepsilon B(T_i)/\\pi``, from its recorded `surface_temperature`.
A facet seen from behind gives `NaN`.

The result is one value per global facet, ready to be placed on an image by a ray-caster such
as `FOVSimulator.generate_image_radiance` (pass ``\\varepsilon = 1`` with the corresponding
[`brightness_temperature`](@ref), or a radiance-taking variant).

# Arguments
- `problem`  : The problem that produced `solution` (shape and emissivities)
- `solution` : Solution with `surface_temperature` recorded (and `roughness_surface_temperature` for rough facets)
- `i_save`   : Index into `solution.output.output_times`
- `d̂`        : Direction from the asteroid to the observer in the body-fixed frame

# Keyword Arguments
- `λ` : Wavelength [m] for the spectral radiance; `nothing` (default) for the total radiance

# Returns
- `L::Vector{Float64}` of length `n_face` [W/m²/sr], or [W/m²/m/sr] when `λ` is given
"""
function directional_radiance(
    problem::SingleAsteroidThermoPhysicalProblem, solution::SingleAsteroidThermoPhysicalSolution,
    i_save::Integer, d̂::StaticVector{3}; λ::Union{Nothing, Real} = nothing,
)
    isnothing(solution.surface_temperature) && throw(ArgumentError(
        "directional_radiance requires save_surface_temperature = true in the output specification"
    ))
    shape = problem.shape
    ε     = problem.thermo_params.emissivity
    d̂     = normalize(d̂)
    T_rough = something(solution.roughness_surface_temperature, Dict{Int, Matrix{Float64}}())

    L = Vector{Float64}(undef, size(solution.surface_temperature, 1))
    for i in eachindex(L)
        if haskey(T_rough, i)
            L[i] = roughness_radiance(shape, i, view(T_rough[i], :, i_save), ε[i], d̂; λ)
        else
            n̂ = shape.face_normals[i]
            L[i] = n̂ ⋅ d̂ > 0 ? ε[i] * _lambert_radiance(solution.surface_temperature[i, i_save], λ) : NaN
        end
    end
    return L
end


"""
    brightness_temperature(problem, solution, i_save, d̂; λ=nothing) -> T_b

Brightness temperature of every global facet towards the observer direction `d̂`: the
temperature of a blackbody whose Lambertian radiance equals the facet's
[`directional_radiance`](@ref), ``B(T_b)/\\pi = L``. Emissivity is not divided out, so a smooth
grey facet at temperature ``T`` has ``T_b = \\varepsilon^{1/4} T`` for the total radiance.

With `λ` the spectral radiance at that wavelength is inverted through the Planck function.
A facet seen from behind gives `NaN`. Returns a vector of length `n_face` [K].
"""
function brightness_temperature(
    problem::SingleAsteroidThermoPhysicalProblem, solution::SingleAsteroidThermoPhysicalSolution,
    i_save::Integer, d̂::StaticVector{3}; λ::Union{Nothing, Real} = nothing,
)
    L = directional_radiance(problem, solution, i_save, d̂; λ)
    return [isnan(l) ? NaN : _brightness_temperature(l, λ) for l in L]
end
