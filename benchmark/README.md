# Performance measurement

Tools for measuring and tracking the performance of the model and its Julia wrappers.
They replace the old `Helios_perf/` package and `src/profile.jl`.

| File | Purpose |
| --- | --- |
| `benchmarks.jl` | The benchmark suite: a BenchmarkTools `SUITE` (PkgBenchmark.jl convention), plus `WORK` (channel-samples per benchmark) for throughput reporting |
| `run.jl` | Runs the suite, prints a table, saves results with provenance metadata, and compares against a saved run |
| `profile.jl` | Statistical profile of one workload with C-level attribution (model code vs. libm vs. GC, by function) |
| `reference.jl` | Output-snapshot check: records all 20 output stages for 13 fixed cases and reports whether later code reproduces them bit-for-bit, to rounding level, or not |
| `results/` | Saved benchmark runs, one directory per run (`results.json` + `meta.toml`); git-ignored |
| `reference/` | Saved output snapshots (about 250 MB each); git-ignored |

## Running

First make sure the C library is current. `deps/build_clib.sh` rebuilds it only when a source file in `deps/model` is newer than the library. A running Julia session uses the rebuilt library on its next model call.

```sh
deps/build_clib.sh            # add --all to also build the 2014 reference library for the tests
```

From a persistent Julia session (for example the Julia MCP server, with the repository root as the environment), which avoids package-load and JIT cost on every run:

```julia
include("benchmark/run.jl")
BenchmarkRunner.main(["--tag", "baseline"])                      # full suite (~1 min)
BenchmarkRunner.main(["--compare", "baseline"])                  # after a change
BenchmarkRunner.main(["--filter", "monaural/", "--quick", "--no-save"])

include("benchmark/profile.jl")
ProfileRunner.main(["--chan", "32", "--dur", "0.5"])
```

The full suite takes longer than the Julia MCP server's default 60-second timeout, so pass `timeout=0` for full runs.

Standalone, from the shell (one-time setup: `julia --project=benchmark -e 'using Pkg; Pkg.instantiate()'`):

```sh
julia --project=benchmark benchmark/run.jl [options]
julia --project=benchmark benchmark/profile.jl [options]
```

### `run.jl` options

| Option | Meaning |
| --- | --- |
| `--filter REGEX` | Only benchmarks whose `group/key` path matches |
| `--tag NAME` | Label stored with the run; usable later as `--compare NAME` |
| `--compare RUN` | Compare against a results directory, a tag (newest run with that tag), or `latest` |
| `--tol X` | Relative tolerance for verdicts (default 0.05) |
| `--quick` | 0.5 s per benchmark; numbers are noisy, so use only as a smoke test |
| `--no-save` | Do not write results |
| `--list` | List benchmark paths |

### Reading the table

- `min` and `median` are wall-clock times per call. `cv` is the coefficient of variation across samples; a high `cv` (more than ~10%) means that row's verdict is not reliable.
- `ns/ch·smp` is the median time divided by the number of simulated channel-samples (ears × channels × samples, including the 20 ms left pad). It is the main size-independent throughput metric. For reference, real time at 100 kHz is 10,000 ns per sample, divided across all channels.
- With `--compare`, `ratio` is the new median divided by the baseline median, and `verdict` applies the `--tol` tolerance to the medians.

## Benchmark groups

| Group | What it measures |
| --- | --- |
| `monaural` | End-to-end `sim_gfc2023` with defaults, scaling in duration and channel count |
| `monaural_nomix` | Same with `moc_width=0.0`; the difference from `monaural` is the O(n_chan²) cross-channel gain-mixing cost |
| `binaural` | End-to-end binaural `sim_gfc2023` |
| `inplace` | `sim_gfc2023!` with a pre-allocated `GFC2023_Mem` and `clean=true` |
| `fractional` | `fractional=true` (includes fractional Gaussian noise synthesis) |
| `powerlaw_exact` | `powerlaw_mode=1` (exact, O(n²) power-law adaptation; used by the regression tests) |
| `components` | `ffGn_native`, memory allocation and zeroing, and fixed per-call overhead |

## Checking that outputs are unchanged

Performance numbers only count if the model still produces the same outputs. `reference.jl` records the outputs of 13 fixed cases once, then compares later code against them:

```julia
include("benchmark/reference.jl")
ReferenceCheck.record()       # once, before optimizing (snapshot "baseline")
ReferenceCheck.compare()      # after each change
```

The cases cover monaural and binaural input; 1, 4, 8 and 32 channels; cross-channel mixing on and off; MOC feedback off; impaired cochlea; cat; fractional noise (seeded); exact power-law adaptation; and the in-place API. Stimuli are generated inside the file from fixed seeds. Recording runs every case twice and fails if the runs differ.

Each case is classified by its worst stage, using the largest difference relative to that stage's peak magnitude:

| Result | Meaning |
| --- | --- |
| `identical` | Bit-for-bit equal. Expected when work is moved or cached without changing any floating-point operation. |
| `rounding` | At most 1e-9. Expected when arithmetic is reordered. For scale, a perturbation of 1e-13 in the MOC gain grows to about 1e-12 through the feedback loop. |
| `small` | At most 1e-6. Inspect before accepting. |
| `CHANGED` | Anything larger, or a change in shape or NaN pattern. The model's behavior changed. |

Binaural cases are insensitive to small changes in the per-ear gain, because the ear-mixing step currently rounds gains to single precision (see "Known issues" in the top-level README).

## Good practice

- Record a tagged baseline before changing anything, and compare against it after each change.
- Run comparisons on the same machine, with the laptop on AC power and other heavy work stopped. `run.jl` warns when a baseline comes from a different host or CPU.
- Results record the commit, whether the tree had uncommitted changes, the Julia version and thread count, the CPU, and the C compiler and flags. Commit before recording a run you want to cite.
- These are performance measurements, not correctness checks. Run the test suite (`include("test/runtests.jl")` in a session) as well.
