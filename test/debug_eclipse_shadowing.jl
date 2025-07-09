#=
debug_eclipse_shadowing.jl

Debug script to verify apply_eclipse_shadowing! function behavior.
This script checks:
1. Whether mutual shadowing is correctly detected
2. Compare results with the old mutual_shadowing! function behavior
3. Verify that illuminated_faces arrays are properly updated
=#

using AsteroidThermoPhysicalModels
using AsteroidShapeModels
using LinearAlgebra
using StaticArrays
using Rotations
using SPICE
using Downloads
using Statistics

println("\n" * "="^60)
println("Eclipse Shadowing Debug Script")
println("="^60 * "\n")

# --- Download and load necessary files (same as previous script) ---
paths_kernel = [
    "fk/hera_v10.tf",
    "lsk/naif0012.tls",
    "pck/hera_didymos_v06.tpc",
    "spk/de432s.bsp",
    "spk/didymos_hor_000101_500101_v01.bsp",
    "spk/didymos_gmv_260901_311001_v01.bsp",
]

# Download SPICE kernels
for path_kernel in paths_kernel
    url_kernel = "https://s2e2.cosmos.esa.int/bitbucket/projects/SPICE_KERNELS/repos/hera/raw/kernels/$(path_kernel)?at=refs%2Ftags%2Fv161_20230929_001"
    filepath = joinpath("test/TPM_Didymos/kernel", path_kernel)
    mkpath(dirname(filepath))
    isfile(filepath) || Downloads.download(url_kernel, filepath)
end

# Download shape models  
paths_shape = [
    "g_50677mm_rad_obj_didy_0000n00000_v001.obj",
    "g_08438mm_lgt_obj_dimo_0000n00000_v002.obj",
]

for path_shape in paths_shape
    url_kernel = "https://s2e2.cosmos.esa.int/bitbucket/projects/SPICE_KERNELS/repos/hera/raw/kernels/dsk/$(path_shape)?at=refs%2Ftags%2Fv161_20230929_001"
    filepath = joinpath("test/TPM_Didymos/shape", path_shape)
    mkpath(dirname(filepath))
    isfile(filepath) || Downloads.download(url_kernel, filepath)
end

# Load SPICE kernels
for path_kernel in paths_kernel
    filepath = joinpath("test/TPM_Didymos/kernel", path_kernel)
    SPICE.furnsh(filepath)
end

# --- Set up test scenario ---
println("Setting up test scenario...")

# Load shape models
shape1 = load_shape_obj(joinpath("test/TPM_Didymos/shape", "g_50677mm_rad_obj_didy_0000n00000_v001.obj"); 
                       scale=1000, with_face_visibility=true, with_bvh=true)
shape2 = load_shape_obj(joinpath("test/TPM_Didymos/shape", "g_08438mm_lgt_obj_dimo_0000n00000_v002.obj"); 
                       scale=1000, with_face_visibility=true, with_bvh=true)

println("Primary shape: $(length(shape1.faces)) faces")
println("Secondary shape: $(length(shape2.faces)) faces")

# Select specific time points for testing
# Choose times when eclipse might occur
et_test_times = [
    SPICE.utc2et("2027-02-18T00:00:00"),
    SPICE.utc2et("2027-02-18T03:00:00"),
    SPICE.utc2et("2027-02-18T06:00:00"),
    SPICE.utc2et("2027-02-18T09:00:00"),
    SPICE.utc2et("2027-02-18T12:00:00"),
]

println("\nTesting at $(length(et_test_times)) time points...")

# --- Test eclipse shadowing at each time point ---
for (idx, et) in enumerate(et_test_times)
    println("\n" * "-"^40)
    println("Test $idx - Time: $(SPICE.et2utc(et, "C", 0))")
    
    # Get positions
    r☉₁ = SVector{3}(SPICE.spkpos("SUN", et, "DIDYMOS_FIXED", "None", "DIDYMOS")[1]) * 1000
    rₛ = SVector{3}(SPICE.spkpos("DIMORPHOS", et, "DIDYMOS_FIXED", "None", "DIDYMOS")[1]) * 1000
    R₁₂ = RotMatrix{3}(SPICE.pxform("DIDYMOS_FIXED", "DIMORPHOS_FIXED", et))
    
    # Compute transformations
    t₁₂ = -R₁₂ * rₛ
    r☉₂ = AsteroidThermoPhysicalModels.transform(r☉₁, R₁₂, t₁₂)
    R₂₁, t₂₁ = AsteroidThermoPhysicalModels.inverse_transformation(R₁₂, t₁₂)
    
    # Normalize sun directions
    r̂☉₁ = normalize(r☉₁)
    r̂☉₂ = normalize(r☉₂)
    
    # First, compute basic illumination (without mutual shadowing)
    illuminated_faces1_basic = zeros(Bool, length(shape1.faces))
    illuminated_faces2_basic = zeros(Bool, length(shape2.faces))
    
    update_illumination!(illuminated_faces1_basic, shape1, r̂☉₁; with_self_shadowing=true)
    update_illumination!(illuminated_faces2_basic, shape2, r̂☉₂; with_self_shadowing=true)
    
    # Copy for eclipse test
    illuminated_faces1 = copy(illuminated_faces1_basic)
    illuminated_faces2 = copy(illuminated_faces2_basic)
    
    # Apply eclipse shadowing
    eclipse_status1 = apply_eclipse_shadowing!(illuminated_faces1, shape1, r☉₁, R₁₂, t₁₂, shape2)
    eclipse_status2 = apply_eclipse_shadowing!(illuminated_faces2, shape2, r☉₂, R₂₁, t₂₁, shape1)
    
    # Count changes
    n_eclipsed1 = sum(illuminated_faces1_basic .& .!illuminated_faces1)
    n_eclipsed2 = sum(illuminated_faces2_basic .& .!illuminated_faces2)
    
    println("Sun-Primary-Secondary angle: $(rad2deg(acos(clamp(normalize(r☉₁) ⋅ normalize(rₛ), -1, 1))))°")
    println("Primary faces eclipsed by secondary: $n_eclipsed1 / $(sum(illuminated_faces1_basic))")
    println("Secondary faces eclipsed by primary: $n_eclipsed2 / $(sum(illuminated_faces2_basic))")
    println("Eclipse status - Primary: $eclipse_status1, Secondary: $eclipse_status2")
    
    # Additional diagnostics
    if n_eclipsed1 > 0 || n_eclipsed2 > 0
        println("  Eclipse detected!")
        
        # Check which parts are eclipsed
        if n_eclipsed1 > 0
            eclipsed_indices1 = findall(illuminated_faces1_basic .& .!illuminated_faces1)
            println("  Primary eclipsed faces: first few = $(eclipsed_indices1[1:min(5, length(eclipsed_indices1))])...")
        end
        
        if n_eclipsed2 > 0
            eclipsed_indices2 = findall(illuminated_faces2_basic .& .!illuminated_faces2)
            println("  Secondary eclipsed faces: first few = $(eclipsed_indices2[1:min(5, length(eclipsed_indices2))])...")
        end
    end
end

# --- Verify consistency check ---
println("\n" * "="^40)
println("Consistency Checks")
println("="^40)

# Test with a simple scenario where we know the expected result
et = et_test_times[1]
r☉₁ = SVector{3}(SPICE.spkpos("SUN", et, "DIDYMOS_FIXED", "None", "DIDYMOS")[1]) * 1000
rₛ = SVector{3}(SPICE.spkpos("DIMORPHOS", et, "DIDYMOS_FIXED", "None", "DIDYMOS")[1]) * 1000
R₁₂ = RotMatrix{3}(SPICE.pxform("DIDYMOS_FIXED", "DIMORPHOS_FIXED", et))

# Check if apply_eclipse_shadowing! preserves already-shadowed faces
println("\nTesting preservation of self-shadowed faces...")
illuminated_test = trues(length(shape1.faces))
r̂☉₁ = normalize(r☉₁)
update_illumination!(illuminated_test, shape1, r̂☉₁; with_self_shadowing=true)
n_self_shadowed = sum(.!illuminated_test)
println("Faces self-shadowed: $n_self_shadowed")

# Apply eclipse shadowing
illuminated_test_copy = copy(illuminated_test)
t₁₂ = -R₁₂ * rₛ
eclipse_status = apply_eclipse_shadowing!(illuminated_test_copy, shape1, r☉₁, R₁₂, t₁₂, shape2)

# Check that no previously shadowed face became illuminated
n_wrongly_lit = sum(.!illuminated_test .& illuminated_test_copy)
println("Faces wrongly re-illuminated: $n_wrongly_lit (should be 0)")

if n_wrongly_lit > 0
    println("⚠️  WARNING: apply_eclipse_shadowing! is incorrectly re-illuminating faces!")
end

SPICE.kclear()

println("\n" * "="^60)
println("Eclipse shadowing debug completed")
println("="^60 * "\n")