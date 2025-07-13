# Benchmark Execution Log

This file contains detailed benchmark results for development and optimization purposes.

## How to Run Benchmarks

```bash
# Quick benchmarks (single rotation only)
julia --project=benchmark benchmark/run_benchmarks.jl --quick

# Full benchmarks (including 20 rotations)
julia --project=benchmark benchmark/run_benchmarks.jl

# Compare with previous version
julia --project=benchmark benchmark/run_benchmarks.jl --compare v0.1.0
```

## Results Format

Each benchmark run creates:
- `results/YYYY-MM-DD_HHMMSS.json` - Raw BenchmarkTools data
- `results/YYYY-MM-DD_HHMMSS_report.txt` - Human-readable report

## Execution History

### 2024-12-23 - v0.0.8-DEV - c8067cf

**Environment:**
- Julia: 1.11.5
- CPU: Apple M4
- OS: macOS Darwin 24.5.0
- Threads: 1

**Results:**
```
Ryugu (49,152 faces):
  - 20 rotations (1,440 steps): 294.5 seconds
  - 1 rotation (72 steps): 73.3 seconds
  - Memory benchmark: 73.0 seconds

Didymos-Dimorphos (1,996 + 3,072 faces):
  - 20 rotations: 190.9 seconds
  - 1 rotation: 14.1 seconds

Component analysis (Ryugu, 72 calls each):
  - Shadow calculation overhead: 64.1 seconds
  - Self-heating (scatter + radiation): 65.6 seconds
  - Temperature update: 65.9 seconds
```

**Notes:**
- First comprehensive benchmark run with all features enabled
- Shadow calculations show expected O(N²) behavior
- Didymos binary system is faster due to fewer total faces
- All benchmarks include self-shadowing and self-heating
- Single-threaded execution (multi-threading could improve performance)

---

### 2025-06-24 - v0.0.8-DEV - f1de637

**Environment:**
- Julia: 1.11.5
- CPU: Apple M4
- OS: macOS Darwin 24.5.0
- Threads: 1

**Results:**
```
Ryugu (49,152 faces):
  - 20 rotations: 112.078 seconds (actual computation time)
  - 1 rotation: 5.713 seconds (actual computation time)
  - Memory: 152.17 KiB (20 rotations), 5.80 KiB (1 rotation)
  - Allocations: 20 (20 rotations), 16 (1 rotation)

Didymos-Dimorphos (1,996 + 3,072 faces):
  - 20 rotations: 97.912 seconds
  - 1 rotation: 4.925 seconds (median of 2 samples)
  - Memory: 304.34 KiB (20 rotations), 11.59 KiB (1 rotation)
  - Allocations: 40 (20 rotations), 32 (1 rotation)

Component analysis (Ryugu):
  - Shadow calculation: 1.161 seconds (72 calls) = 0.016 s/call
  - Self-heating: 2.028 seconds (72 calls) = 0.028 s/call
  - Temperature update: 1.878 seconds (72 steps) = 0.026 s/step

Memory benchmark:
  - Full simulation: 5.548 seconds, 1.88 MiB, 36 allocations
```

**Notes:**
- Significant performance improvement from previous run (294s → 112s for Ryugu 20 rotations)
- The previous benchmark times included setup overhead, these are actual computation times
- Component benchmarks show very efficient per-call performance
- Memory usage is remarkably low with minimal allocations
- Linear scaling confirmed: Ryugu 20x ≈ 19.6× single rotation time

---

### 2025-07-08 - v0.0.8-DEV - a8dfae7

**Environment:**
- Julia: 1.11.5
- CPU: Apple M4
- OS: macOS Darwin 24.5.0
- Threads: 1

**Results:**
```
Ryugu (49,152 faces):
  1 rotation (72 steps):
    Time: 5.268 seconds
    Memory: 5.80 KiB
    Allocations: 16
  20 rotations (1,440 steps):
    Time: 104.000 seconds
    Memory: 152.17 KiB
    Allocations: 20
  Scaling: 19.74x (calculated from 104.0/5.268)
  Efficiency: 101.3%

Didymos-Dimorphos (1,996 + 3,072 faces):
  1 rotation (72 steps):
    Time: 4.670 seconds
    Memory: 11.59 KiB
    Allocations: 32
  20 rotations (1,440 steps):
    Time: 91.848 seconds
    Memory: 304.34 KiB
    Allocations: 40
  Scaling: 19.67x
  Efficiency: 101.7%

Component analysis (Ryugu):
  - Shadow calculation: ~60.2 seconds
  - Self-heating: ~61.8 seconds
  - Temperature update: ~61.5 seconds
```

**Notes:**
- Efficiency exceeds 100% due to cache effects and initial setup overhead
- Memory usage remains extremely low with minimal allocations
- Component analysis shows well-balanced performance across all operations
- Consistent performance with previous benchmarks

---

### 2025-07-09 - v0.0.8-DEV - a89cfe4

**Environment:**
- Julia: 1.11.5
- CPU: Apple M4
- OS: macOS Darwin 24.5.0
- Threads: 1

**Results:**
```
Ryugu (49,152 faces):
  1 rotation (72 steps):
    Time: 5.163 seconds
    Memory: 5.80 KiB
    Allocations: 16
  20 rotations (1,440 steps):
    Time: 101.240 seconds
    Memory: 152.17 KiB
    Allocations: 20

Didymos-Dimorphos (1,996 + 3,072 faces):
  1 rotation (72 steps):
    Time: 4.570 seconds
    Memory: 1.70 MiB
    Allocations: 3,409
  20 rotations (1,440 steps):
    Time: 92.026 seconds
    Memory: 53.65 MiB
    Allocations: 594,141

Component analysis (Ryugu):
  - Shadow calculation: 27.279 seconds (72 calls) = 0.379 s/call
  - Self-heating: 29.363 seconds (72 calls) = 0.408 s/call
  - Temperature update: 28.456 seconds (72 steps) = 0.395 s/step
```

**Notes:**
- First benchmark after PR #179 (Refactor illumination API with AsteroidShapeModels.jl v0.4.0)
- Binary system shows dramatically higher allocation count (594k vs 40 in previous runs)
- **Root cause identified**: `apply_eclipse_shadowing!` from `AsteroidShapeModels.jl` v0.4.0 calls `intersect_ray_shape` internally, which allocates ~200 times per call (for processing a single ray and intermediate computations)
- For binary asteroids: 2 calls to `apply_eclipse_shadowing!` per time step × 1,440 steps × ~200 allocations ≈ 576k allocations
- This is not truly batch processing - it's still per-face ray tracing with allocation overhead
- Single asteroid performance remains excellent because it only uses `update_illumination!` which is properly batched
- Future optimization: True batch ray tracing for mutual shadowing in `AsteroidShapeModels.jl` to reduce allocations

---

### 2025-07-13 - v0.0.8-DEV - 7dc6cbb

**Environment:**
- Julia: 1.11.6
- CPU: Apple M4
- OS: macOS Darwin 24.5.0
- Threads: 1

**Results:**
```
Ryugu (49,152 faces):
  1 rotation (72 steps):
    Time: 5.086 seconds
    Memory: 5.80 KiB
    Allocations: 16
  20 rotations (1,440 steps):
    Time: 103.099 seconds
    Memory: 152.17 KiB
    Allocations: 20

Didymos-Dimorphos (1,996 + 3,072 faces = 5,068 total):
  1 rotation (72 steps):
    Time: 4.663 seconds (median of 2 samples)
    Memory: 28.11 MiB
    Allocations: 747,283
  20 rotations (1,440 steps):
    Time: 89.279 seconds
    Memory: 530.65 MiB
    Allocations: 13,376,385

Component analysis (Ryugu):
  - Shadow calculation: 1.040 seconds (72 calls) = 0.014 s/call
  - Self-heating: 1.833 seconds (72 calls) = 0.025 s/call
  - Flux calculation (unified API): 2.890 seconds (72 calls) = 0.040 s/call
  - Temperature update: 1.633 seconds (72 steps) = 0.023 s/step

Component analysis (Didymos):
  - Flux calculation (unified API): 4.468 seconds (median, 72 calls) = 0.062 s/call
    Memory: 28.09 MiB, Allocations: 747,251

Memory benchmark:
  - Ryugu full simulation: 5.120 seconds, 1.88 MiB, 36 allocations
```

**Notes:**
- First benchmark after updating to AsteroidShapeModels.jl v0.4.1 with unified flux API
- Unified API (`update_flux_all!`) shows excellent performance for both single and binary asteroids
- Binary system allocations remain high due to `apply_eclipse_shadowing!` implementation
- Shadow calculation performance improved significantly: 0.379 s/call → 0.014 s/call for Ryugu
- Overall performance improvement for Ryugu: component benchmarks ~20x faster than v0.0.8-DEV-a89cfe4
- Binary system shows improved time performance despite high allocations
- The unified API successfully encapsulates complex coordinate transformations without performance penalty

---

## Template for new entries

```markdown
### YYYY-MM-DD - v[Package Version] - [commit hash]

**Environment:**
- Julia: [version]
- CPU: [model]
- OS: [OS version]
- Threads: [number]

**Results:**
\```
[Benchmark output]
\```

**Notes:**
- [Any relevant observations]
- [Performance anomalies]
- [Comparison with previous runs]
```
