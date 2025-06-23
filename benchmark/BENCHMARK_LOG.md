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

### [Date] - [Git commit hash]

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