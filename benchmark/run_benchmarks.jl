#!/usr/bin/env julia

# Script to run benchmarks for AsteroidThermoPhysicalModels.jl

using Pkg
Pkg.activate(@__DIR__)

# Ensure required packages are available
required_packages = [
    "AsteroidShapeModels",
    "AsteroidThermoPhysicalModels",
    "BenchmarkTools",
    "Dates",
    "Downloads",
    "LinearAlgebra",
    "Rotations",
    "SPICE",
    "StaticArrays",
    "Statistics"
]

println("Checking required packages...")
for pkg in required_packages
    if !haskey(Pkg.project().dependencies, pkg)
        println("Adding $pkg...")
        Pkg.add(pkg)
    end
end

# Include benchmark definitions
include("benchmarks.jl")

# Import Dates for timestamp formatting
using Dates

# Run benchmarks
println("\n" * "="^60)
println("  Running AsteroidThermoPhysicalModels.jl Benchmarks")
println("="^60)

# Download required data first
println("\nDownloading benchmark data...")
download_benchmark_data()

# Run selected benchmarks
println("\nRunning benchmarks...")
println("This may take several minutes for the full 20-rotation simulations.")

# You can customize which benchmarks to run
# For quick test, use single rotation benchmarks:
# suite_to_run = BenchmarkGroup()
# suite_to_run["ryugu"] = BenchmarkGroup()
# suite_to_run["ryugu"]["single_rotation"] = SUITE["ryugu"]["single_rotation"]
# results = run(suite_to_run, verbose=true)

# For full benchmarks:
results = run(SUITE, verbose=true)

# Save results
timestamp = Dates.format(Dates.now(), "yyyy-mm-dd_HHMMSS")
results_file = joinpath(@__DIR__, "results", "benchmark_results_$(timestamp).json")
mkpath(dirname(results_file))
BenchmarkTools.save(results_file, results)

println("\nResults saved to: $results_file")

# Print summary
print_summary(results)

# Get package version
try
    pkg_version = Pkg.project().version
    if isnothing(pkg_version)
        pkg_version = "dev"
    end
catch
    pkg_version = "unknown"
end

# Export human-readable report
report_file = joinpath(@__DIR__, "results", "benchmark_report_$(timestamp).txt")
open(report_file, "w") do io
    # Redirect output to file
    original_stdout = stdout
    redirect_stdout(io)
    
    println("AsteroidThermoPhysicalModels.jl Benchmark Report")
    println("Version: v$(pkg_version)")
    println("Generated: $(Dates.now())")
    println("="^60)
    
    for (group_name, group) in results
        println("\n## $group_name")
        println("-"^40)
        for (bench_name, bench) in group
            println("\n### $bench_name")
            show(io, MIME("text/plain"), bench)
            println()
        end
    end
    
    redirect_stdout(original_stdout)
end

println("Report saved to: $report_file")
println("\nBenchmark complete!")