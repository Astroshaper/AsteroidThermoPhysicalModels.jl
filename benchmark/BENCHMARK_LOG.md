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

### [Date] - v0.0.8-DEV - [Git commit hash]

**Environment:**
- Julia: [version]
- CPU: [model]
- OS: [OS version]
- Threads: [number]

**Results:**
```
[Paste benchmark output here]
```

**Notes:**
- [Any relevant observations]
- [Performance anomalies]
- [Comparison with previous runs]

---

### Template for new entries

```markdown
### YYYY-MM-DD - v[Package Version] - [commit hash]

**Environment:**
- Julia: 
- CPU: 
- OS: 
- Threads: 

**Results:**
\```
[Benchmark output]
\```

**Notes:**
- 
```