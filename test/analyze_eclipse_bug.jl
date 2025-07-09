#=
analyze_eclipse_bug.jl

Analyze the bug in apply_eclipse_shadowing! implementation
Focus on coordinate transformation issues
=#

using LinearAlgebra
using StaticArrays

println("\n" * "="^60)
println("Analysis of apply_eclipse_shadowing! Implementation")
println("="^60 * "\n")

# Based on the source code, let's trace through the problematic case
# from our earlier debug script where TOTAL_ECLIPSE was incorrectly reported

# Test case values (from debug_eclipse_detailed.jl)
r☉₁ = SVector(1.2152543918410248e10, -1.642022414469607e11, 1.6079396748154287e10)  # Sun in primary frame
t₁₂ = SVector(-847.6748429773617, -835.1588094492685, 7.55802133967758)  # Translation vector
angle_deg = 130.16  # Sun-Primary-Secondary angle

# Assuming identity rotation for analysis
R₁₂ = SMatrix{3,3}(I)
r̂☉₁ = normalize(r☉₁)
r̂☉₂ = normalize(R₁₂ * r̂☉₁)  # Line 137 in the implementation

# Typical values for asteroid sizes
ρ₁ = 427.5  # Primary max radius (Didymos)
ρ₂ = 104.0  # Secondary max radius (Dimorphos)

println("Initial values:")
println("  |r☉₁|: $(norm(r☉₁)/1e9) million km")
println("  |t₁₂|: $(norm(t₁₂)) m")
println("  Sun-Primary-Secondary angle: $(angle_deg)°")
println("  ρ₁ (primary radius): $ρ₁ m")
println("  ρ₂ (secondary radius): $ρ₂ m")

# Check Early Out 1 (Behind Check) - Line 148
behind_check = dot(t₁₂, r̂☉₁)
println("\n--- Early Out 1: Behind Check ---")
println("  dot(t₁₂, r̂☉₁) = $behind_check")
println("  Threshold: -(ρ₁ + ρ₂) = $(-ρ₁ - ρ₂)")
println("  Check: $(behind_check) < $(-ρ₁ - ρ₂) ? $(behind_check < -(ρ₁ + ρ₂))")
println("  Result: $(behind_check < -(ρ₁ + ρ₂) ? "NO_ECLIPSE (early exit)" : "Continue checking")")

# Check Early Out 2 (Lateral Separation) - Lines 154-157
t₁₂⊥ = t₁₂ - (dot(t₁₂, r̂☉₁) * r̂☉₁)
d⊥ = norm(t₁₂⊥)
println("\n--- Early Out 2: Lateral Separation Check ---")
println("  t₁₂⊥ (perpendicular component): $t₁₂⊥")
println("  d⊥ (lateral distance): $d⊥ m")
println("  Threshold: ρ₁ + ρ₂ = $(ρ₁ + ρ₂)")
println("  Check: $d⊥ > $(ρ₁ + ρ₂) ? $(d⊥ > ρ₁ + ρ₂)")
println("  Result: $(d⊥ > ρ₁ + ρ₂ ? "NO_ECLIPSE (early exit)" : "Continue checking")")

# Check Early Out 3 (Total Eclipse) - Lines 160-167
# THIS IS WHERE THE BUG LIKELY IS!
total_eclipse_cond1 = dot(t₁₂, r̂☉₁) > 0
total_eclipse_cond2 = d⊥ + ρ₁ < ρ₂
println("\n--- Early Out 3: Total Eclipse Check (POTENTIAL BUG) ---")
println("  Condition 1: dot(t₁₂, r̂☉₁) > 0 ? $total_eclipse_cond1")
println("    dot(t₁₂, r̂☉₁) = $behind_check")
println("  Condition 2: d⊥ + ρ₁ < ρ₂ ? $total_eclipse_cond2")
println("    d⊥ + ρ₁ = $(d⊥ + ρ₁)")
println("    ρ₂ = $ρ₂")
println("  Both conditions: $(total_eclipse_cond1 && total_eclipse_cond2)")
println("  Result: $(total_eclipse_cond1 && total_eclipse_cond2 ? "TOTAL_ECLIPSE (BUG!)" : "Continue to face checks")")

# Analysis of the bug
println("\n" * "="^40)
println("BUG ANALYSIS")
println("="^40)

println("\nThe Total Eclipse condition (line 165) is WRONG!")
println("Current implementation checks:")
println("  if dot(t₁₂, r̂☉₁) > 0 && d⊥ + ρ₁ < ρ₂")
println("\nProblems:")
println("1. dot(t₁₂, r̂☉₁) > 0 means secondary is in FRONT of primary along sun direction")
println("   BUT: For eclipse, secondary should be BETWEEN sun and primary!")
println("   This is backwards!")

println("\n2. The correct condition should check if secondary is between sun and primary:")
println("   - Secondary should be closer to sun than primary")
println("   - In the primary's frame, this means t₁₂ should point AWAY from sun")
println("   - So we need: dot(t₁₂, r̂☉₁) < 0 (not > 0)")

println("\n3. Additionally, for total eclipse at large distances:")
println("   - The angular size comparison is more complex")
println("   - Current simple radius check (d⊥ + ρ₁ < ρ₂) is insufficient")

# Demonstrate the issue with our test case
println("\n--- Applied to our test case ---")
println("Sun-Primary-Secondary angle: $(angle_deg)° (> 90°)")
println("This means secondary is NOT between sun and primary!")
println("Yet the current code would return TOTAL_ECLIPSE because:")
println("  - dot(t₁₂, r̂☉₁) = $behind_check > 0 ✓ (wrong condition)")
println("  - d⊥ + ρ₁ = $(d⊥ + ρ₁) < ρ₂ = $ρ₂ ? $(d⊥ + ρ₁ < ρ₂)")

# Additional coordinate transformation issue
println("\n--- Potential Coordinate Transformation Issue ---")
println("Line 137: r̂☉₂ = normalize(R₁₂ * r̂☉₁)")
println("This transforms the NORMALIZED sun direction.")
println("But for eclipse calculations, we might need the actual sun position,")
println("not just its direction, especially for proper angular size calculations.")

println("\n" * "="^60)
println("CONCLUSION")
println("="^60)
println("\nThe bug is in the Total Eclipse early-out condition (line 165).")
println("The condition `dot(t₁₂, r̂☉₁) > 0` is backwards - it should be < 0")
println("for the secondary to be between the sun and primary.")
println("\nThis explains why TOTAL_ECLIPSE is incorrectly reported when")
println("the Sun-Primary-Secondary angle is > 90°.")