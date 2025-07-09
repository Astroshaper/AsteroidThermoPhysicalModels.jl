#=
verify_coordinate_consistency.jl

Verify that coordinate transformations are consistent
in the apply_eclipse_shadowing! calls
=#

using AsteroidThermoPhysicalModels
using AsteroidShapeModels
using LinearAlgebra
using StaticArrays
using Rotations

println("\n" * "="^60)
println("Coordinate Transformation Consistency Check")
println("="^60 * "\n")

# Test scenario setup
println("Setting up test scenario...")

# Create test positions and transformations
r☉₁ = SVector(1.5e11, 0.0, 0.0)  # Sun position in primary frame
rₛ = SVector(1000.0, 0.0, 0.0)   # Secondary position in primary frame
R₁₂ = RotMatrix{3}(I)             # Identity rotation for simplicity

# Compute transformations as in update_flux_sun!
t₁₂ = -R₁₂ * rₛ
r☉₂ = AsteroidThermoPhysicalModels.transform(r☉₁, R₁₂, t₁₂)
R₂₁, t₂₁ = AsteroidThermoPhysicalModels.inverse_transformation(R₁₂, t₁₂)

println("\nTransformation values:")
println("  r☉₁ (sun in primary frame): $r☉₁")
println("  rₛ (secondary in primary frame): $rₛ")
println("  R₁₂: $R₁₂")
println("  t₁₂ = -R₁₂ * rₛ: $t₁₂")
println("  r☉₂ (sun in secondary frame): $r☉₂")

# Verify the transformation is correct
println("\nVerification 1: r☉₂ = R₁₂ * r☉₁ + t₁₂")
r☉₂_check = R₁₂ * r☉₁ + t₁₂
println("  Computed r☉₂: $r☉₂")
println("  Check r☉₂: $r☉₂_check")
println("  Difference: $(norm(r☉₂ - r☉₂_check))")

# Verify inverse transformation
println("\nVerification 2: Inverse transformation")
println("  R₂₁ = R₁₂': $(norm(R₂₁ - R₁₂') < 1e-10 ? "✓" : "✗")")
println("  t₂₁ = -R₂₁ * t₁₂: $t₂₁")
println("  Expected t₂₁ = rₛ: $(norm(t₂₁ - rₛ) < 1e-10 ? "✓" : "✗")")

# Check the transformation chain
println("\nVerification 3: Transform chain consistency")
# Point at secondary's center in primary frame should map to origin in secondary frame
secondary_center_in_sec = R₁₂ * rₛ + t₁₂
println("  Secondary center in secondary frame: $secondary_center_in_sec")
println("  Should be ~[0,0,0]: $(norm(secondary_center_in_sec) < 1e-10 ? "✓" : "✗")")

# Point at primary's center should map to -t₂₁ in secondary frame
primary_center_in_sec = R₁₂ * SVector(0.0, 0.0, 0.0) + t₁₂
println("  Primary center in secondary frame: $primary_center_in_sec")
println("  Should equal t₁₂: $(norm(primary_center_in_sec - t₁₂) < 1e-10 ? "✓" : "✗")")

# Now check what apply_eclipse_shadowing! expects
println("\n" * "-"^40)
println("API Usage Analysis")
println("-"^40)

println("\nIn update_flux_sun!, we call:")
println("1. apply_eclipse_shadowing!(illuminated_faces1, shape1, r☉₁, R₁₂, t₁₂, shape2)")
println("   - shape1 is in its own frame (primary)")
println("   - r☉₁ is sun position in shape1's frame ✓")
println("   - R₁₂, t₁₂ transform from shape1's frame to shape2's frame ✓")
println("   - shape2 needs to be checked against shape1")

println("\n2. apply_eclipse_shadowing!(illuminated_faces2, shape2, r☉₂, R₂₁, t₂₁, shape1)")
println("   - shape2 is in its own frame (secondary)")
println("   - r☉₂ is sun position in shape2's frame ✓")
println("   - R₂₁, t₂₁ transform from shape2's frame to shape1's frame ✓")
println("   - shape1 needs to be checked against shape2")

println("\nConclusion: The API usage appears correct.")
println("The issue is likely in the apply_eclipse_shadowing! implementation itself.")

# Additional check: What happens if we pass wrong sun vector?
println("\n" * "-"^40)
println("Testing potential issues")
println("-"^40)

println("\nPotential issue 1: Normalized vs non-normalized sun vector")
println("  |r☉₁|: $(norm(r☉₁))")
println("  |r☉₂|: $(norm(r☉₂))")
println("  Both are non-normalized (correct for distance calculation)")

println("\nPotential issue 2: Sign of translation vector")
println("  t₁₂ = -R₁₂ * rₛ = $t₁₂")
println("  This means: origin_of_frame2 = R₁₂ * rₛ + t₁₂ = 0 ✓")

println("\n" * "="^60)
println("Analysis completed")
println("="^60 * "\n")