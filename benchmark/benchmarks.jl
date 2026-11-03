# Benchmark suite for caplab_binaural_efferent_model.
#
# Defines `SUITE::BenchmarkGroup` (the PkgBenchmark.jl / AirspeedVelocity.jl convention) plus
# `WORK`, which records how many channel-samples of model work each benchmark performs so
# that results can be reported as a size-independent throughput (ns per channel-sample).
#
# Run it with benchmark/run.jl (see benchmark/README.md) rather than including this file
# directly. Stimuli are seeded, so every run benchmarks exactly the same inputs.

using BenchmarkTools
using Random
using caplab_binaural_efferent_model
const BEM = caplab_binaural_efferent_model

const FS = 100e3
const DUR_PAD_LEFT = 0.02   # wrapper default; adds 2000 samples to every simulation
const SUITE = BenchmarkGroup()
const WORK = Dict{String,Float64}()   # "group/key" => channel-samples (ears x chans x samples)

# ~60 dB SPL Gaussian noise (0.02 Pa rms); seeded so all runs see identical input
stim(dur; seed=1) = 0.02 .* randn(Xoshiro(seed), round(Int, dur * FS))
stim2(dur) = [stim(dur; seed=1), stim(dur; seed=2)]
cfs(n) = n == 1 ? [1000.0] : collect(exp.(range(log(125.0), log(16e3), n)))
nsamp(dur; pad=DUR_PAD_LEFT) = round(Int, dur * FS) + floor(Int, pad * FS)

"""
    add!(group, key, bench; work=nothing, seconds=2.0)

Register `bench` as `SUITE[group][key]`, with a time budget of `seconds`. `work` is the
number of channel-samples simulated, used to report ns/channel-sample.
"""
function add!(group::String, key::String, bench::BenchmarkTools.Benchmark; work=nothing, seconds=2.0)
    haskey(SUITE, group) || (SUITE[group] = BenchmarkGroup())
    bench.params.seconds = seconds
    bench.params.evals = 1          # every call is ≥ 1 ms; never batch evaluations
    bench.params.samples = 10_000   # bounded by `seconds`
    bench.params.gctrial = true
    SUITE[group][key] = bench
    work === nothing || (WORK["$group/$key"] = work)
    return bench
end

label(n, dur) = "chan=$(lpad(n, 3)) dur=$(dur)s"

# ------------------------------------------------------------------------------------------
# End-to-end monaural simulations (default parameters, incl. cross-channel MOC mixing).
# The main throughput benchmarks: scaling in duration (rows 1-2) and channels (rows 3-5).
# ------------------------------------------------------------------------------------------
for (n, dur, sec) in [(1, 0.1, 2.0), (1, 1.0, 3.0), (8, 0.1, 2.0), (32, 0.1, 3.0), (128, 0.02, 5.0)]
    x, c = stim(dur), cfs(n)
    add!("monaural", label(n, dur), @benchmarkable(sim_gfc2023($x, $c));
         work=n * nsamp(dur), seconds=sec)
end

# Same, without cross-channel MOC mixing (moc_width=0): the difference from "monaural"
# isolates the O(n_chan^2) per-sample gain-mixing cost.
for (n, dur, sec) in [(32, 0.1, 3.0), (128, 0.02, 5.0)]
    x, c = stim(dur), cfs(n)
    add!("monaural_nomix", label(n, dur), @benchmarkable(sim_gfc2023($x, $c; moc_width=0.0));
         work=n * nsamp(dur), seconds=sec)
end

# ------------------------------------------------------------------------------------------
# End-to-end binaural simulations
# ------------------------------------------------------------------------------------------
for (n, dur, sec) in [(1, 0.1, 2.0), (8, 0.1, 2.0), (32, 0.1, 3.0)]
    x, c = stim2(dur), cfs(n)
    add!("binaural", label(n, dur), @benchmarkable(sim_gfc2023($x, $c));
         work=2 * n * nsamp(dur), seconds=sec)
end

# ------------------------------------------------------------------------------------------
# In-place API with pre-allocated memory (compare against "monaural" of the same size to
# see what pre-allocation buys). clean=true is required for correct results on reuse.
# ------------------------------------------------------------------------------------------
for (n, dur, sec) in [(8, 0.1, 2.0), (32, 0.1, 3.0)]
    x, c = stim(dur), cfs(n)
    mem = GFC2023_Mem(x, c)
    add!("inplace", label(n, dur), @benchmarkable(sim_gfc2023!($mem, $x, $c; clean=true));
         work=n * nsamp(dur), seconds=sec)
end
for (n, dur, sec) in [(8, 0.1, 2.0), (32, 0.1, 3.0)]
    x, c = stim2(dur), cfs(n)
    mem = GFC2023_Mem(x, c)
    add!("inplace_binaural", label(n, dur), @benchmarkable(sim_gfc2023!($mem, $x, $c; clean=true));
         work=2 * n * nsamp(dur), seconds=sec)
end

# ------------------------------------------------------------------------------------------
# Fractional Gaussian noise (fractional=true adds ffGn synthesis for 2 fiber types/channel)
# ------------------------------------------------------------------------------------------
for (n, dur, sec) in [(1, 1.0, 3.0), (8, 0.1, 2.0)]
    x, c = stim(dur), cfs(n)
    add!("fractional", label(n, dur), @benchmarkable(sim_gfc2023($x, $c; fractional=true));
         work=n * nsamp(dur), seconds=sec)
end

# ------------------------------------------------------------------------------------------
# Exact power-law adaptation (powerlaw_mode=1), O(n^2) in duration; kept short. Not used by
# default but used by the regression test suite.
# ------------------------------------------------------------------------------------------
let dur = 0.05
    x = stim(dur)
    add!("powerlaw_exact", label(1, dur), @benchmarkable(sim_gfc2023($x, [1000.0]; powerlaw_mode=1));
         work=nsamp(dur), seconds=3.0)
end

# ------------------------------------------------------------------------------------------
# Components and fixed per-call overheads
# ------------------------------------------------------------------------------------------
add!("components", "ffGn_native 1s", @benchmarkable(BEM.ffGn_native($(nsamp(1.0)), 1 / FS, 0.9, 1.0, 100.0)))
add!("components", "GFC2023_Mem chan=32 dur=0.1s", @benchmarkable(GFC2023_Mem($(stim(0.1)), $(cfs(32)))))
let mem = GFC2023_Mem(stim(0.1), cfs(32))
    add!("components", "zero_state! chan=32 dur=0.1s", @benchmarkable(zero_state!($mem)))
end
# Tiny simulation without padding: dominated by fixed per-call costs (keyword defaults,
# allocation, library lookup, C-side setup) rather than per-sample work
let x = stim(0.001)
    add!("components", "call overhead (100 samples)", @benchmarkable(sim_gfc2023($x, [1000.0]; dur_pad_left=0.0));
         work=length(x))
end
