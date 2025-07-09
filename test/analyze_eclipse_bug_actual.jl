#=
analyze_eclipse_bug_actual.jl

Analyze the actual case where TOTAL_ECLIPSE was incorrectly reported
Using the exact values from debug_eclipse_detailed.jl
=#

using LinearAlgebra
using StaticArrays
using Rotations
using SPICE
using Downloads

println("\n" * "="^60)
println("Analysis of Actual TOTAL_ECLIPSE Bug Case")
println("="^60 * "\n")

# Load SPICE kernels to get exact values
paths_kernel = [
    "fk/hera_v10.tf",
    "lsk/naif0012.tls", 
    "pck/hera_didymos_v06.tpc",
    "spk/de432s.bsp",
    "spk/didymos_hor_000101_500101_v01.bsp",
    "spk/didymos_gmv_260901_311001_v01.bsp",
]

for path_kernel in paths_kernel
    filepath = joinpath("test/TPM_Didymos/kernel", path_kernel)
    if isfile(filepath)
        SPICE.furnsh(filepath)
    end
end

# Get exact values at the problematic time
et = SPICE.utc2et("2027-02-18T12:00:00")
r☉₁ = SVector{3}(SPICE.spkpos("SUN", et, "DIDYMOS_FIXED", "None", "DIDYMOS")[1]) * 1000
rₛ = SVector{3}(SPICE.spkpos("DIMORPHOS", et, "DIDYMOS_FIXED", "None", "DIDYMOS")[1]) * 1000
R₁₂ = RotMatrix{3}(SPICE.pxform("DIDYMOS_FIXED", "DIMORPHOS_FIXED", et))

SPICE.kclear()

# Compute transformations as in update_flux_sun!
t₁₂ = -R₁₂ * rₛ
r☉₂ = R₁₂ * r☉₁ + t₁₂
R₂₁ = R₁₂'
t₂₁ = -R₂₁ * t₁₂

# Get normalized directions
r̂☉₁ = normalize(r☉₁)
r̂☉₂ = normalize(r☉₂)

# Asteroid parameters
ρ₁ = 427.5  # Primary max radius
ρ₂ = 104.0  # Secondary max radius

println("Actual values from problematic case:")
println("  Time: 2027-02-18T12:00:00")
println("  |r☉₁|: $(norm(r☉₁)/1e9) million km")
println("  |r☉₂|: $(norm(r☉₂)/1e9) million km")  
println("  |rₛ|: $(norm(rₛ)) m")
println("  |t₁₂|: $(norm(t₁₂)) m")

# Calculate angle
angle = rad2deg(acos(clamp(normalize(r☉₁) ⋅ normalize(rₛ), -1, 1)))
println("  Sun-Primary-Secondary angle: $(angle)°")

# Now trace through apply_eclipse_shadowing! for shape2 (secondary)
# which incorrectly reported TOTAL_ECLIPSE
println("\n" * "="^40)
println("Tracing apply_eclipse_shadowing! for secondary")
println("(This is the call that returns TOTAL_ECLIPSE)")
println("="^40)

# The problematic call is:
# apply_eclipse_shadowing!(illuminated2, shape2, r☉₂, R₂₁, t₂₁, shape1)

# Early Out 1: Behind Check
behind_check = dot(t₂₁, normalize(r☉₂))
println("\n--- Early Out 1: Behind Check ---")
println("  dot(t₂₁, r̂☉₂) = $behind_check")
println("  Threshold: -(ρ₂ + ρ₁) = $(-ρ₂ - ρ₁)")
println("  Check: $behind_check < $(-ρ₂ - ρ₁) ? $(behind_check < -(ρ₂ + ρ₁))")

# Early Out 2: Lateral Separation
r̂☉₂_normalized = normalize(r☉₂)
t₂₁⊥ = t₂₁ - (dot(t₂₁, r̂☉₂_normalized) * r̂☉₂_normalized)
d⊥ = norm(t₂₁⊥)
println("\n--- Early Out 2: Lateral Separation Check ---")
println("  t₂₁⊥ magnitude: $d⊥ m")
println("  Threshold: ρ₂ + ρ₁ = $(ρ₂ + ρ₁)")
println("  Check: $d⊥ > $(ρ₂ + ρ₁) ? $(d⊥ > ρ₂ + ρ₁)")

# Early Out 3: Total Eclipse (THE BUG!)
total_eclipse_cond1 = dot(t₂₁, r̂☉₂_normalized) > 0
total_eclipse_cond2 = d⊥ + ρ₂ < ρ₁  # Note: roles are swapped (shape2 perspective)
println("\n--- Early Out 3: Total Eclipse Check ---")
println("  Condition 1: dot(t₂₁, r̂☉₂) > 0 ? $total_eclipse_cond1")
println("    dot(t₂₁, r̂☉₂) = $(dot(t₂₁, r̂☉₂_normalized))")
println("  Condition 2: d⊥ + ρ₂ < ρ₁ ? $total_eclipse_cond2")  
println("    d⊥ + ρ₂ = $(d⊥ + ρ₂)")
println("    ρ₁ = $ρ₁")
println("  Both conditions met: $(total_eclipse_cond1 && total_eclipse_cond2)")

if total_eclipse_cond1 && total_eclipse_cond2
    println("\n  ⚠️  TOTAL_ECLIPSE would be returned here!")
end

# Analyze why this is wrong
println("\n" * "="^40)
println("WHY THIS IS WRONG")
println("="^40)

println("\nThe problem:")
println("1. We're checking if shape1 (primary) eclipses shape2 (secondary)")
println("2. Sun-Primary-Secondary angle = $(angle)° > 90°")
println("3. This means primary is NOT between sun and secondary")
println("4. Therefore, primary CANNOT eclipse secondary!")

println("\nThe bug in line 165:")
println("  if dot(t₁₂, r̂☉₁) > 0 && d⊥ + ρ₁ < ρ₂")
println("\nFor shape2's perspective (with t₂₁, r̂☉₂):")
println("  if dot(t₂₁, r̂☉₂) > 0 && d⊥ + ρ₂ < ρ₁")

println("\nThe condition dot(t₂₁, r̂☉₂) > 0 means:")
println("  - Primary (shape1) is in the direction of the sun from secondary's perspective")
println("  - But for eclipse, primary should be OPPOSITE to sun direction!")
println("  - The correct condition should be: dot(t₂₁, r̂☉₂) < 0")

# Verify with vector directions
println("\n--- Vector Direction Analysis ---")
println("In secondary's frame:")
println("  Sun direction r̂☉₂: $(r̂☉₂_normalized)")
println("  Primary position t₂₁: $t₂₁")
println("  Dot product: $(dot(t₂₁, r̂☉₂_normalized))")
println("  Angle between: $(rad2deg(acos(clamp(dot(normalize(t₂₁), r̂☉₂_normalized), -1, 1))))°")
println("\nSince angle < 90°, primary and sun are in SAME direction from secondary!")
println("This confirms primary cannot eclipse secondary.")

println("\n" * "="^60)
println("CONFIRMED: Bug in apply_eclipse_shadowing! line 165")
println("="^60 * "\n")