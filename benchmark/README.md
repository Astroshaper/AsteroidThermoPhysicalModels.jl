# AsteroidThermoPhysicalModels.jl Benchmarks

This directory contains performance benchmarks for tracking the computational efficiency of the package across different versions.

## Running Benchmarks

### Quick Start

```bash
julia benchmark/run_benchmarks.jl
```

### Using PkgBenchmark.jl

For comparing performance between different commits or branches:

```julia
using PkgBenchmark
using AsteroidThermoPhysicalModels

# Run benchmarks on current state
results = benchmarkpkg(AsteroidThermoPhysicalModels)

# Compare with main branch
judge_results = judge(AsteroidThermoPhysicalModels, "main")

# Compare with specific commit
judge_results = judge(AsteroidThermoPhysicalModels, "abc123")
```

## Benchmark Suite

The benchmark suite includes:

### 1. **Ryugu Simulations**
- Single rotation (baseline)
- 20 rotations (full benchmark)
- ~49,000 faces shape model
- Self-shadowing and self-heating enabled

### 2. **Didymos-Dimorphos Binary System**
- Single rotation (baseline)
- 20 rotations (full benchmark)
- Primary: ~50,000 faces
- Secondary: ~8,000 faces
- Mutual shadowing and heating enabled

### 3. **Component Analysis**
- Shadow calculation overhead
- Self-heating calculations
- Temperature update performance
- Memory allocation patterns

## Expected Performance

On a modern desktop CPU (e.g., Intel i7 or AMD Ryzen 7):

- **Ryugu (49k faces, 1 rotation)**: ~30-60 seconds
- **Ryugu (49k faces, 20 rotations)**: ~10-20 minutes
- **Didymos Binary (1 rotation)**: ~45-90 seconds
- **Didymos Binary (20 rotations)**: ~15-30 minutes

## Benchmark Data

The benchmarks automatically download required SPICE kernels and shape models:
- Ryugu shape model from JAXA/Hayabusa2
- Didymos system models from ESA/Hera
- Ephemeris data for realistic sun positions

## Results

Benchmark results are saved in:
- `benchmark/results/benchmark_results_YYYY-MM-DD_HHMMSS.json` - Raw data
- `benchmark/results/benchmark_report_YYYY-MM-DD_HHMMSS.txt` - Human-readable report

## Performance Tracking

To track performance over time:

1. Run benchmarks before making changes
2. Make your optimizations
3. Run benchmarks again
4. Compare results using the provided comparison functions

## Optimization Targets

Key areas for optimization based on profiling:
1. **Shadow calculations** (isilluminated function) - O(N²) complexity
2. **Visibility graph queries** - Frequent allocations
3. **Heat conduction solvers** - Lack of SIMD optimizations
4. **Memory allocations** - Unnecessary array copies