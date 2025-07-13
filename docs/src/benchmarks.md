# Performance and Benchmarks

This page provides performance expectations and optimization guidelines for AsteroidThermoPhysicalModels.jl.

## Expected Performance

The following benchmarks were performed on Apple M4 (macOS, single-threaded):

### Single Asteroid (Ryugu)
- **Shape complexity**: 49,152 faces
- **1 rotation** (72 time steps): ~5.1 seconds
- **20 rotations** (1,440 time steps): ~103 seconds (1.7 minutes)
- **Memory usage**: 152 KiB (20 rotations), 5.8 KiB (1 rotation)
- **Allocations**: 20 (20 rotations), 16 (1 rotation)
- **With shadows and self-heating enabled**

### Binary System (Didymos-Dimorphos)
- **Primary**: 1,996 faces, **Secondary**: 3,072 faces (5,068 total)
- **1 rotation** (72 time steps): ~4.7 seconds (based on Dimorphos' rotation period)
- **20 rotations** (1,440 time steps): ~89 seconds (1.5 minutes)
- **Memory usage**: 531 MiB (20 rotations), 28.1 MiB (1 rotation)
- **Allocations**: 13.4M (20 rotations), 747k (1 rotation)
- **With mutual shadowing and heating enabled**
- **Note**: Rotation counts refer to Dimorphos (secondary) rotation periods

### Component Performance (per time step)
- **Shadow calculations**: ~0.014 seconds (1.0s for 72 steps) - **27x faster with v0.4.1**
- **Self-heating**: ~0.025 seconds (1.8s for 72 steps)
- **Temperature update**: ~0.023 seconds (1.6s for 72 steps)
- **Unified flux calculation**: ~0.040 seconds (2.9s for 72 steps) for Ryugu
- **Unified flux calculation**: ~0.062 seconds (4.5s for 72 steps) for Didymos

*Note: Component benchmarks show total time for 72 calls in isolation*

*Note: Performance may vary depending on CPU architecture. Intel/AMD processors may show different characteristics.*

## Performance Considerations

### Computational Complexity

The main computational bottlenecks are:

1. **Shadow calculations**: O(N²) where N is the number of faces
   - Dominates computation time for large shape models
   - Can be disabled with `SELF_SHADOWING = false` for faster computation

2. **Self-heating**: O(N×M) where M is the average number of visible faces
   - Uses precomputed visibility graph
   - Can be disabled with `SELF_HEATING = false`

3. **Temperature solver**: O(N×D) where D is the number of depth layers
   - Typically not a bottleneck unless using many depth layers

### Memory Usage

Memory scales approximately as:
- Shape model: O(N)
- Temperature array: O(N×D)
- Visibility graph: O(N×M) where M is average visible faces per face

For Ryugu (49k faces, 41 depth layers):
- Expected memory usage: 2-3 GB
- Peak during visibility computation: 4-5 GB

## Optimization Tips

### 1. Disable Features for Faster Computation

```julia
# Fastest configuration (no shadows or self-heating)
stpm = SingleAsteroidTPM(shape, thermo_params;
    SELF_SHADOWING = false,
    SELF_HEATING   = false
)
```

### 2. Reduce Shape Complexity

For preliminary calculations, use a simplified shape model.

### 3. Adjust Time Resolution

Larger time steps can be used for initial calculations:
```julia
# Use fewer steps per rotation
n_step_in_cycle = 36  # Instead of 72
```

## Version Performance History

### v0.1.0 (Current)
- Updated to `AsteroidShapeModels.jl` v0.4.1 with critical bug fixes
- Implemented unified flux update API (`update_flux_all!`)
- Fixed coordinate transformation bug in binary systems
- **Performance**: Ryugu 20 rotations in ~103s, Didymos in ~89s
- **Memory**: Ryugu uses 152 KiB with minimal allocations (20)
- **Allocations**: Binary system shows high allocation count (13.4M) due to `apply_eclipse_shadowing!` implementation
- **Shadow calculations**: 27x faster (0.379 → 0.014 s/call)
- **Latest benchmark**: 2025-07-13 (Apple M4)

### v0.0.7
- Memory optimizations in flux calculations
- Basic shadow and self-heating calculations

### Earlier versions
- Initial implementation
- Basic thermophysical modeling

## Benchmarking Your Configuration

To benchmark your specific use case:

```julia
using BenchmarkTools
using AsteroidThermoPhysicalModels

# Your setup
stpm = # ... your model
ephem = # ... your ephemerides

# Benchmark
@benchmark run_TPM!($stpm, $ephem, Float64[], Int[]; show_progress=false)
```

For detailed benchmarking, see the `benchmark/` directory in the source repository.
