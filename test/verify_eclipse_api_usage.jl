#=
verify_eclipse_api_usage.jl

Verify the correct usage of apply_eclipse_shadowing! function
by examining the function signature and our usage in update_flux_sun!
=#

using Pkg
using AsteroidShapeModels

println("\n" * "="^60)
println("Verifying apply_eclipse_shadowing! API usage")
println("="^60 * "\n")

# Check the method signatures for apply_eclipse_shadowing!
println("Available methods for apply_eclipse_shadowing!:")
methods_list = methods(apply_eclipse_shadowing!)
for m in methods_list
    println("\n", m)
end

# Let's also check what AsteroidShapeModels exports
println("\n\nExported names from AsteroidShapeModels:")
exported_names = names(AsteroidShapeModels)
eclipse_related = filter(n -> occursin("eclipse", lowercase(string(n))), exported_names)
println("Eclipse-related exports: ", eclipse_related)

# Check for documentation
println("\n\nTrying to access documentation...")
try
    # This might not work outside REPL, but let's try
    doc = @doc apply_eclipse_shadowing!
    println(doc)
catch e
    println("Could not access documentation: ", e)
end

# Check the package version
println("\n\nAsteroidShapeModels.jl version:")
pkg_info = Pkg.dependencies()
for (uuid, info) in pkg_info
    if info.name == "AsteroidShapeModels"
        println("  Version: ", info.version)
        println("  Source: ", info.source)
        break
    end
end

println("\n" * "="^60)
println("Verification completed")
println("="^60 * "\n")