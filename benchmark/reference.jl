# Output-snapshot check: record the model's outputs for a fixed set of cases, then compare
# later versions of the code against that record. Used to verify that optimizations do not
# change model behavior. It checks consistency with the recorded version, not scientific
# correctness.
#
#   include("benchmark/reference.jl")            # from a persistent session (e.g. Julia MCP)
#   ReferenceCheck.record()                      # record snapshot "baseline"
#   ReferenceCheck.compare()                     # compare current code against "baseline"
#   ReferenceCheck.compare("baseline"; cases=["mono_8ch"])   # subset of cases
#
#   julia --project=benchmark benchmark/reference.jl record [NAME]    # standalone
#   julia --project=benchmark benchmark/reference.jl compare [NAME]
#
# Snapshots are stored in benchmark/reference/NAME/ (git-ignored): one serialized file of
# outputs per case plus meta.toml recording the commit and environment they came from.
# Serialized files are specific to the Julia version that wrote them.
#
# For every case and output stage, compare() reports the largest absolute difference and
# the largest difference relative to the stage's peak magnitude, and classifies the case:
#   identical  bit-for-bit equal (expected when work is only moved, cached, or reordered
#              without changing any floating-point operation)
#   rounding   relative difference <= 1e-9 (expected when arithmetic is reordered)
#   small      relative difference <= 1e-6 (inspect: plausible for reordered arithmetic
#              inside the feedback loops, but worth a look)
#   CHANGED    anything larger, or a change in shape or NaN pattern

module ReferenceCheck

using Dates
using Printf
using Random
using Serialization
using TOML
using caplab_binaural_efferent_model

const ROOT = normpath(joinpath(@__DIR__, ".."))
const REF_DIR = joinpath(@__DIR__, "reference")
const FS = 100e3

const STAGES = ["control", "c1", "c2", "ihc", "expon_hsr", "sout1_hsr", "sout2_hsr", "syn_hsr",
                "expon_lsr", "sout1_lsr", "sout2_lsr", "syn_lsr", "hsr", "lsr", "cn", "ic",
                "mocwdr", "mocic", "gain", "gainpostmix"]

# ------------------------------------------------------------------------------------------
# Stimuli (self-contained so that changes in other packages cannot alter the inputs)
# ------------------------------------------------------------------------------------------
rms_pa(level_db) = 20e-6 * 10^(level_db / 20)

function ramp!(x; dur_ramp=0.01)
    n = round(Int, dur_ramp * FS)
    for i in 1:n
        w = 0.5 * (1 - cos(pi * (i - 1) / n))
        x[i] *= w
        x[end-i+1] *= w
    end
    return x
end

tone(f, level, dur) = ramp!(sqrt(2) * rms_pa(level) .* sin.(2pi * f .* (0:round(Int, dur * FS)-1) ./ FS))
noise(level, dur; seed) = ramp!(rms_pa(level) .* randn(Xoshiro(seed), round(Int, dur * FS)))
cfs(n) = n == 1 ? [1000.0] : collect(exp.(range(log(125.0), log(16e3), n)))

run_inplace(x, c) = sim_gfc2023!(GFC2023_Mem(x, c), x, c; clean=true)

# ------------------------------------------------------------------------------------------
# Cases: name => zero-argument function returning the raw sim_gfc2023 output
# ------------------------------------------------------------------------------------------
const CASES = Pair{String,Function}[
    "mono_1ch_tone" => () -> sim_gfc2023(tone(1000.0, 60.0, 0.2), [1000.0]),
    "mono_1ch_loud" => () -> sim_gfc2023(noise(85.0, 0.3; seed=11), [4000.0]),
    "mono_8ch" => () -> sim_gfc2023(noise(70.0, 0.2; seed=12), cfs(8)),
    "mono_8ch_nomix" => () -> sim_gfc2023(noise(70.0, 0.2; seed=12), cfs(8); moc_width=0.0),
    "mono_8ch_nomoc" => () -> sim_gfc2023(noise(70.0, 0.2; seed=12), cfs(8); moc_weight=0.0),
    "mono_32ch" => () -> sim_gfc2023(noise(70.0, 0.1; seed=13), cfs(32)),
    "mono_4ch_impaired" => () -> sim_gfc2023(noise(75.0, 0.2; seed=14), cfs(4);
                                             cohc=fill(0.5, 4), cihc=fill(0.8, 4)),
    "mono_4ch_cat" => () -> sim_gfc2023(noise(70.0, 0.2; seed=15), cfs(4); species="cat"),
    "mono_4ch_fractional" => () -> (Random.seed!(16); sim_gfc2023(noise(70.0, 0.2; seed=16), cfs(4); fractional=true)),
    "mono_1ch_powerlaw_exact" => () -> sim_gfc2023(tone(1000.0, 60.0, 0.05), [1000.0]; powerlaw_mode=1),
    "mono_8ch_inplace" => () -> run_inplace(noise(70.0, 0.2; seed=12), cfs(8)),
    "binaural_8ch" => () -> sim_gfc2023([noise(70.0, 0.2; seed=17), noise(65.0, 0.2; seed=18)], cfs(8)),
    "binaural_1ch_tone" => () -> sim_gfc2023([tone(1000.0, 70.0, 0.2), tone(1000.0, 50.0, 0.2)], [1000.0]),
]

# Convert raw output into Dict(stage => [ear matrix (channel x time), ...])
function to_snapshot(out)
    snap = Dict{String,Vector{Matrix{Float64}}}()
    for (k, stage) in enumerate(STAGES)
        s = out[k]
        ears = s[1] isa AbstractVector{<:AbstractVector} ? s : [s]   # binaural: [ear][chan][time]
        snap[stage] = [permutedims(reduce(hcat, ear)) for ear in ears]
    end
    return snap
end

# ------------------------------------------------------------------------------------------
# Recording
# ------------------------------------------------------------------------------------------
gitcmd(args...) = try readchomp(Cmd(`git -C $ROOT $args`)) catch; "unknown" end

"""
    record(name="baseline"; overwrite=false)

Run every case twice (failing if the two runs differ, which would indicate nondeterminism)
and save the outputs as snapshot `name`.
"""
function record(name::AbstractString="baseline"; overwrite::Bool=false)
    dir = joinpath(REF_DIR, name)
    if isdir(dir)
        overwrite || error("snapshot \"$name\" already exists at $dir; pass overwrite=true to replace it")
        rm(dir; recursive=true)
    end
    mkpath(dir)
    t0 = time()
    for (case, f) in CASES
        a, b = to_snapshot(f()), to_snapshot(f())
        for stage in STAGES
            isequal(a[stage], b[stage]) || error("case $case is nondeterministic: stage $stage differs between two identical runs")
        end
        serialize(joinpath(dir, case * ".jls"), a)
    end
    meta = Dict{String,Any}(
        "name" => name,
        "timestamp" => string(now()),
        "git_commit" => gitcmd("rev-parse", "HEAD"),
        "git_dirty" => !isempty(gitcmd("status", "--porcelain", "--untracked-files=no")),
        "julia_version" => string(VERSION),
        "cases" => first.(CASES),
    )
    open(io -> TOML.print(io, meta), joinpath(dir, "meta.toml"), "w")
    size_mb = sum(filesize(joinpath(dir, f)) for f in readdir(dir)) / 1024^2
    println("Recorded snapshot \"$name\": $(length(CASES)) cases, $(round(size_mb; digits=1)) MB, ",
            "$(round(time() - t0; digits=1)) s, commit $(meta["git_commit"][1:min(8, end)])",
            meta["git_dirty"] ? " (dirty)" : "")
    return dir
end

# ------------------------------------------------------------------------------------------
# Comparison
# ------------------------------------------------------------------------------------------
const CLASSES = ["identical", "rounding", "small", "CHANGED"]

function classify(ref::Vector{Matrix{Float64}}, new::Vector{Matrix{Float64}})
    if length(ref) != length(new) || any(size(r) != size(n) for (r, n) in zip(ref, new))
        return (class="CHANGED", abs=NaN, rel=NaN, note="shape changed")
    end
    isequal(ref, new) && return (class="identical", abs=0.0, rel=0.0, note="")
    max_abs, max_rel, nan_mismatch = 0.0, 0.0, false
    for (r, n) in zip(ref, new)
        nan_mismatch |= any(isnan.(r) .!= isnan.(n))
        d = maximum(abs.(ifelse.(isnan.(r) .| isnan.(n), 0.0, r .- n)); init=0.0)
        scale = maximum(abs.(filter(!isnan, r)); init=0.0)
        max_abs = max(max_abs, d)
        max_rel = max(max_rel, scale > 0 ? d / scale : (d > 0 ? Inf : 0.0))
    end
    nan_mismatch && return (class="CHANGED", abs=max_abs, rel=max_rel, note="NaN pattern changed")
    class = max_rel <= 1e-9 ? "rounding" : max_rel <= 1e-6 ? "small" : "CHANGED"
    return (class=class, abs=max_abs, rel=max_rel, note="")
end

"""
    compare(name="baseline"; cases=nothing, verbose=false)

Run the cases (all, or the names in `cases`) with the current code and compare against
snapshot `name`. Prints one line per case and, for non-identical cases (or every case if
`verbose`), one line per affected stage. Returns the worst class found.
"""
function compare(name::AbstractString="baseline"; cases=nothing, verbose::Bool=false)
    dir = joinpath(REF_DIR, name)
    isdir(dir) || error("no snapshot \"$name\" in $REF_DIR; record one with ReferenceCheck.record(\"$name\")")
    meta = TOML.parsefile(joinpath(dir, "meta.toml"))
    meta["julia_version"] == string(VERSION) ||
        @warn "snapshot was recorded with Julia $(meta["julia_version"]); serialized files may not load"
    selected = cases === nothing ? CASES : filter(p -> first(p) in cases, CASES)
    println("Comparing against snapshot \"$name\" (commit $(meta["git_commit"][1:min(8, end)])",
            meta["git_dirty"] ? ", dirty" : "", ", recorded $(meta["timestamp"][1:16]))")
    println(rpad("case", 26), rpad("result", 11), lpad("max rel", 10), lpad("max abs", 11), "  worst stage")
    println("-"^72)
    worst = 1
    for (case, f) in selected
        file = joinpath(dir, case * ".jls")
        if !isfile(file)
            println(rpad(case, 26), "(not in snapshot)")
            continue
        end
        ref, new = deserialize(file), to_snapshot(f())
        per_stage = [(stage, classify(ref[stage], new[stage])) for stage in STAGES]
        ci = maximum(findfirst(==(c.class), CLASSES) for (_, c) in per_stage)
        worst = max(worst, ci)
        stage_worst, cw = argmax(((s, c),) -> (findfirst(==(c.class), CLASSES), c.rel), per_stage)
        println(rpad(case, 26), rpad(CLASSES[ci], 11), lpad(@sprintf("%.1e", cw.rel), 10),
                lpad(@sprintf("%.1e", cw.abs), 11), "  ", ci == 1 ? "-" : stage_worst)
        if verbose || ci > 1
            for (s, c) in per_stage
                (verbose || c.class != "identical") || continue
                println("    ", rpad(s, 22), rpad(c.class, 11), lpad(@sprintf("%.1e", c.rel), 10),
                        lpad(@sprintf("%.1e", c.abs), 11), isempty(c.note) ? "" : "  " * c.note)
            end
        end
    end
    println("\nOverall: ", CLASSES[worst])
    return CLASSES[worst]
end

end # module ReferenceCheck

if abspath(PROGRAM_FILE) == @__FILE__
    cmd = isempty(ARGS) ? "compare" : ARGS[1]
    name = length(ARGS) >= 2 ? ARGS[2] : "baseline"
    if cmd == "record"
        ReferenceCheck.record(name)
    elseif cmd == "compare"
        exit(ReferenceCheck.compare(name) == "CHANGED" ? 1 : 0)
    else
        error("usage: reference.jl record|compare [NAME]")
    end
end
