#=
debug_eclipse_detailed.jl

Detailed investigation of the eclipse shadowing behavior.
Focus on the case where secondary shows TOTAL_ECLIPSE.
=#

using AsteroidThermoPhysicalModels
using AsteroidShapeModels
using LinearAlgebra
using StaticArrays
using Rotations
using SPICE
using Downloads

println("\n" * "="^60)
println("Detailed Eclipse Investigation")
println("="^60 * "\n")

# Load kernels properly
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

# Load shapes
shape1 = load_shape_obj(joinpath("test/TPM_Didymos/shape", "g_50677mm_rad_obj_didy_0000n00000_v001.obj"); 
                       scale=1000, with_face_visibility=true, with_bvh=true)
shape2 = load_shape_obj(joinpath("test/TPM_Didymos/shape", "g_08438mm_lgt_obj_dimo_0000n00000_v002.obj"); 
                       scale=1000, with_face_visibility=true, with_bvh=true)

# Focus on the problematic time (Test 5 from previous script)
et = SPICE.utc2et("2027-02-18T12:00:00")
println("Investigating time: $(SPICE.et2utc(et, "C", 0))")

# Get positions and transformations
r☉₁ = SVector{3}(SPICE.spkpos("SUN", et, "DIDYMOS_FIXED", "None", "DIDYMOS")[1]) * 1000
rₛ = SVector{3}(SPICE.spkpos("DIMORPHOS", et, "DIDYMOS_FIXED", "None", "DIDYMOS")[1]) * 1000
R₁₂ = RotMatrix{3}(SPICE.pxform("DIDYMOS_FIXED", "DIMORPHOS_FIXED", et))

t₁₂ = -R₁₂ * rₛ
r☉₂ = AsteroidThermoPhysicalModels.transform(r☉₁, R₁₂, t₁₂)
R₂₁, t₂₁ = AsteroidThermoPhysicalModels.inverse_transformation(R₁₂, t₁₂)

println("\nPositions and distances:")
println("  Sun distance from primary: $(norm(r☉₁)/1e9) million km")
println("  Sun distance from secondary: $(norm(r☉₂)/1e9) million km") 
println("  Primary-secondary distance: $(norm(rₛ)/1000) km")
println("  Sun-Primary-Secondary angle: $(rad2deg(acos(clamp(normalize(r☉₁) ⋅ normalize(rₛ), -1, 1))))°")

# Check basic illumination first
r̂☉₁ = normalize(r☉₁)
r̂☉₂ = normalize(r☉₂)

illuminated_faces2_basic = zeros(Bool, length(shape2.faces))
update_illumination!(illuminated_faces2_basic, shape2, r̂☉₂; with_self_shadowing=true)
println("\nSecondary basic illumination: $(sum(illuminated_faces2_basic)) faces lit (out of $(length(shape2.faces)))")

# Now test eclipse shadowing step by step
println("\nTesting eclipse shadowing on secondary...")

# Method 1: Direct call as in the code
illuminated_faces2_test1 = copy(illuminated_faces2_basic)
eclipse_status1 = apply_eclipse_shadowing!(illuminated_faces2_test1, shape2, r☉₂, R₂₁, t₂₁, shape1)
n_eclipsed1 = sum(illuminated_faces2_basic .& .!illuminated_faces2_test1)
println("Method 1 - Faces eclipsed: $n_eclipsed1, Status: $eclipse_status1")

# Check if the issue is with the sun position passed to apply_eclipse_shadowing!
println("\nDiagnostic checks:")
println("  r☉₂ (sun in secondary frame): $r☉₂")
println("  |r☉₂|: $(norm(r☉₂))")
println("  R₂₁ * r☉₂ + t₂₁ (should equal r☉₁): $(R₂₁ * r☉₂ + t₂₁)")
println("  Verification error: $(norm(R₂₁ * r☉₂ + t₂₁ - r☉₁))")

# Test with normalized vs non-normalized sun vector
println("\nTesting with different sun vector inputs...")

# Test 2: With normalized sun vector (incorrect usage?)
illuminated_faces2_test2 = copy(illuminated_faces2_basic)
eclipse_status2 = apply_eclipse_shadowing!(illuminated_faces2_test2, shape2, normalize(r☉₂), R₂₁, t₂₁, shape1)
n_eclipsed2 = sum(illuminated_faces2_basic .& .!illuminated_faces2_test2)
println("With normalized sun: Faces eclipsed: $n_eclipsed2, Status: $eclipse_status2")

# Check shape bounds
println("\nShape bounds check:")
println("  Primary max radius: $(maximum_radius(shape1)) m")
println("  Secondary max radius: $(maximum_radius(shape2)) m")

# Manual eclipse check
# When Sun-Primary-Secondary angle is > 90°, secondary cannot be in primary's shadow
angle = rad2deg(acos(clamp(normalize(r☉₁) ⋅ normalize(rₛ), -1, 1)))
if angle > 90
    println("\n⚠️  Sun-Primary-Secondary angle is $(angle)° > 90°")
    println("   Secondary CANNOT be in primary's shadow!")
    println("   But apply_eclipse_shadowing! returns TOTAL_ECLIPSE")
    println("   This indicates a BUG in the eclipse detection algorithm")
end

# Additional check: Primary eclipse on secondary
# The primary should never cause total eclipse when angle > 90°
println("\nVerifying eclipse geometry...")
# Vector from secondary to primary (in secondary frame)
primary_pos_in_sec = R₂₁ * rₛ + t₂₁  # Should be -t₂₁ since secondary is at origin
println("  Primary position in secondary frame: $primary_pos_in_sec")
println("  Expected (should be -t₂₁): $(-t₂₁)")
println("  Difference: $(norm(primary_pos_in_sec - (-t₂₁)))")

SPICE.kclear()

println("\n" * "="^60)
println("Investigation completed")
println("="^60 * "\n")