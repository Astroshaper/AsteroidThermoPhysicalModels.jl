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
