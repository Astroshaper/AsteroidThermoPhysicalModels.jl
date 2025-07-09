#=
debug_coordinate_transformation.jl

Debug script to verify coordinate transformation equivalence between PR #178/#179 (before) and current implementation (after).
This script checks that:
1. ephem.sun2 (old) == R₁₂ * r☉₁ + t₁₂ (new)
2. ephem.S2P (old) == R₁₂' (new)
=#

using AsteroidThermoPhysicalModels
using LinearAlgebra
using StaticArrays
using Rotations
using SPICE
using Downloads
using Statistics

println("\n" * "="^60)
println("Coordinate Transformation Debug Script")
println("="^60 * "\n")

# --- Download and load SPICE kernels (same as TPM_Didymos test) ---
paths_kernel = [
    "fk/hera_v10.tf",
    "lsk/naif0012.tls",
    "pck/hera_didymos_v06.tpc",
    "spk/de432s.bsp",
    "spk/didymos_hor_000101_500101_v01.bsp",
    "spk/didymos_gmv_260901_311001_v01.bsp",
]

for path_kernel in paths_kernel
    url_kernel = "https://s2e2.cosmos.esa.int/bitbucket/projects/SPICE_KERNELS/repos/hera/raw/kernels/$(path_kernel)?at=refs%2Ftags%2Fv161_20230929_001"
    filepath = joinpath("test/TPM_Didymos/kernel", path_kernel)
    mkpath(dirname(filepath))
    isfile(filepath) || Downloads.download(url_kernel, filepath)
end

for path_kernel in paths_kernel
    filepath = joinpath("test/TPM_Didymos/kernel", path_kernel)
    SPICE.furnsh(filepath)
end

# --- Set up ephemerides ---
P₁ = SPICE.convrt(2.2593, "hours", "seconds")  # Rotation period of Didymos
P₂ = SPICE.convrt(11.93 , "hours", "seconds")  # Rotation period of Dimorphos

n_cycle = 1  # Just one cycle for debugging
n_step_in_cycle = 10  # Fewer steps for quick debugging

et_begin = SPICE.utc2et("2027-02-18T00:00:00")
et_end   = et_begin + P₂ * n_cycle
et_range = range(et_begin, et_end; length=n_step_in_cycle*n_cycle+1)

# --- Compute ephemerides using the OLD method (simulating pre-PR behavior) ---
println("Computing ephemerides using OLD method (pre-PR #179)...")

ephem_old = (
    time = collect(et_range),
    sun1 = [SVector{3}(SPICE.spkpos("SUN", et, "DIDYMOS_FIXED", "None", "DIDYMOS")[1]) * 1000 for et in et_range],
    sun2 = [SVector{3}(SPICE.spkpos("SUN", et, "DIMORPHOS_FIXED", "None", "DIMORPHOS")[1]) * 1000 for et in et_range],
    sec  = [SVector{3}(SPICE.spkpos("DIMORPHOS", et, "DIDYMOS_FIXED", "None", "DIDYMOS")[1]) * 1000 for et in et_range],
    P2S  = [RotMatrix{3}(SPICE.pxform("DIDYMOS_FIXED", "DIMORPHOS_FIXED", et)) for et in et_range],
    S2P  = [RotMatrix{3}(SPICE.pxform("DIMORPHOS_FIXED", "DIDYMOS_FIXED", et)) for et in et_range],
)

# --- Compute ephemerides using the NEW method (current implementation) ---
println("Computing ephemerides using NEW method (post-PR #179)...")

ephem_new = (
    time = collect(et_range),
    sun  = [SVector{3}(SPICE.spkpos("SUN", et, "DIDYMOS_FIXED", "None", "DIDYMOS")[1]) * 1000 for et in et_range],
    sec  = [SVector{3}(SPICE.spkpos("DIMORPHOS", et, "DIDYMOS_FIXED", "None", "DIDYMOS")[1]) * 1000 for et in et_range],
    P2S  = [RotMatrix{3}(SPICE.pxform("DIDYMOS_FIXED", "DIMORPHOS_FIXED", et)) for et in et_range],
)

SPICE.kclear()

# --- Verification 1: Check that sun1 (old) == sun (new) ---
println("\n--- Verification 1: sun1 (old) vs sun (new) ---")
sun_diff = [norm(ephem_old.sun1[i] - ephem_new.sun[i]) for i in eachindex(et_range)]
println("Maximum difference: $(maximum(sun_diff)) m")
println("Mean difference: $(mean(sun_diff)) m")
@assert maximum(sun_diff) < 1e-10 "sun1 and sun should be identical!"

# --- Verification 2: Check S2P (old) == P2S' (new) ---
println("\n--- Verification 2: S2P (old) vs P2S' (new) ---")
s2p_diff = Float64[]
for i in eachindex(et_range)
    R₁₂ = ephem_new.P2S[i]
    R₂₁_old = ephem_old.S2P[i]
    R₂₁_new = R₁₂'
    
    # Compare rotation matrices element by element
    diff = norm(R₂₁_old - R₂₁_new)
    push!(s2p_diff, diff)
end
println("Maximum difference in rotation matrices: $(maximum(s2p_diff))")
println("Mean difference in rotation matrices: $(mean(s2p_diff))")
@assert maximum(s2p_diff) < 1e-10 "S2P should equal P2S'!"

# --- Verification 3: Check sun2 (old) == R₁₂ * sun + t₁₂ (new) ---
println("\n--- Verification 3: sun2 (old) vs transform(sun, R₁₂, t₁₂) (new) ---")
sun2_diff = Float64[]
sun2_computed = Vector{SVector{3, Float64}}()

for i in eachindex(et_range)
    # Old method: directly from SPICE
    r☉₂_old = ephem_old.sun2[i]
    
    # New method: compute using transformation
    r☉₁ = ephem_new.sun[i]
    rₛ = ephem_new.sec[i]
    R₁₂ = ephem_new.P2S[i]
    
    # Translation vector calculation (as in current code)
    t₁₂ = -R₁₂ * rₛ
    
    # Transform sun position to secondary frame
    r☉₂_new = AsteroidThermoPhysicalModels.transform(r☉₁, R₁₂, t₁₂)
    push!(sun2_computed, r☉₂_new)
    
    # Compare
    diff = norm(r☉₂_old - r☉₂_new)
    push!(sun2_diff, diff)
    
    if i <= 3  # Print first few for inspection
        println("\nTime step $i:")
        println("  sun2 (old):      $(r☉₂_old)")
        println("  sun2 (computed): $(r☉₂_new)")
        println("  Difference:      $diff m")
    end
end

println("\n--- Summary of sun2 differences ---")
println("Maximum difference: $(maximum(sun2_diff)) m")
println("Mean difference: $(mean(sun2_diff)) m")
println("Median difference: $(median(sun2_diff)) m")

# Check if differences are within acceptable tolerance
tolerance = 1e-6  # 1 micrometer tolerance
if maximum(sun2_diff) > tolerance
    println("\n⚠️  WARNING: Sun position differences exceed tolerance!")
    println("This could explain the unnatural temperature distribution.")
else
    println("\n✓ Sun positions match within tolerance.")
end

# --- Additional check: Verify the transformation formula ---
println("\n--- Verification 4: Mathematical consistency check ---")
for i in [1, div(length(et_range), 2), length(et_range)]  # Check beginning, middle, end
    r☉₁ = ephem_new.sun[i]
    rₛ = ephem_new.sec[i]
    R₁₂ = ephem_new.P2S[i]
    t₁₂ = -R₁₂ * rₛ
    
    # Method 1: Using transform function
    r☉₂_method1 = AsteroidThermoPhysicalModels.transform(r☉₁, R₁₂, t₁₂)
    
    # Method 2: Manual calculation
    r☉₂_method2 = R₁₂ * r☉₁ + t₁₂
    
    # Method 3: Alternative formulation
    r☉₂_method3 = R₁₂ * (r☉₁ - rₛ)
    
    println("\nTime step $i - Transformation methods comparison:")
    println("  Method 1 (transform function): $(r☉₂_method1)")
    println("  Method 2 (R*r + t):           $(r☉₂_method2)")
    println("  Method 3 (R*(r-rs)):          $(r☉₂_method3)")
    println("  Diff 1-2: $(norm(r☉₂_method1 - r☉₂_method2))")
    println("  Diff 1-3: $(norm(r☉₂_method1 - r☉₂_method3))")
    println("  Diff 2-3: $(norm(r☉₂_method2 - r☉₂_method3))")
end

# --- Check inverse transformation consistency ---
println("\n--- Verification 5: Inverse transformation consistency ---")
for i in [1, length(et_range)]
    R₁₂ = ephem_new.P2S[i]
    rₛ = ephem_new.sec[i]
    t₁₂ = -R₁₂ * rₛ
    
    # Compute inverse transformation
    R₂₁, t₂₁ = AsteroidThermoPhysicalModels.inverse_transformation(R₁₂, t₁₂)
    
    # Check if R₂₁ * R₁₂ = I
    identity_check = R₂₁ * R₁₂
    identity_error = norm(identity_check - I)
    
    # Check if we can recover rₛ
    rₛ_in_sec_frame = R₁₂ * rₛ + t₁₂  # Should be zero vector
    rₛ_recovered = R₂₁ * rₛ_in_sec_frame + t₂₁  # Should equal rₛ
    
    println("\nTime step $i:")
    println("  Identity check error: $identity_error")
    println("  rₛ in secondary frame: $(norm(rₛ_in_sec_frame)) (should be ~0)")
    println("  rₛ recovery error: $(norm(rₛ_recovered - rₛ)) m")
end

println("\n" * "="^60)
println("Debug script completed")
println("="^60 * "\n")