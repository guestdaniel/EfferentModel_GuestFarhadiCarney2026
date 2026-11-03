# Run the benchmark suite, save results with provenance metadata, and optionally compare
# against a previous run. See benchmark/README.md.
#
#   julia --project=benchmark benchmark/run.jl [options]      # standalone, from the shell
#
#   include("benchmark/run.jl")                               # from a persistent session
#   BenchmarkRunner.main(["--compare", "baseline"])           # (e.g. the Julia MCP server)
#
# Options:
#   --filter REGEX   Only run benchmarks whose "group/key" path matches REGEX
#   --tag NAME       Label stored with the run (e.g. baseline); usable with --compare
#   --compare RUN    Compare against RUN: a results directory, a tag, or "latest"
#   --tol X          Relative tolerance for regression/improvement verdicts (default 0.05)
#   --quick          Short time budget per benchmark (smoke test; numbers are noisy)
#   --no-save        Do not write results to benchmark/results
#   --list           List benchmark paths and exit

module BenchmarkRunner

using BenchmarkTools
using Dates
using Printf
using Statistics
using TOML

const BENCH_DIR = @__DIR__
const ROOT = normpath(joinpath(BENCH_DIR, ".."))
const RESULTS_DIR = joinpath(BENCH_DIR, "results")

include(joinpath(BENCH_DIR, "benchmarks.jl"))

# ------------------------------------------------------------------------------------------
# Argument parsing
# ------------------------------------------------------------------------------------------
function parse_args(args)
    opts = Dict{String,Any}("filter" => nothing, "tag" => "", "compare" => nothing, "tol" => 0.05,
                            "quick" => false, "save" => true, "list" => false)
    i = 1
    while i <= length(args)
        a = args[i]
        if a in ("--filter", "--tag", "--compare", "--tol")
            i < length(args) || error("missing value for $a")
            opts[a[3:end]] = a == "--tol" ? parse(Float64, args[i+1]) : args[i+1]
            i += 2
        elseif a == "--quick"
            opts["quick"] = true; i += 1
        elseif a == "--no-save"
            opts["save"] = false; i += 1
        elseif a == "--list"
            opts["list"] = true; i += 1
        else
            error("unknown argument: $a (see header of benchmark/run.jl)")
        end
    end
    return opts
end

# ------------------------------------------------------------------------------------------
# Provenance metadata
# ------------------------------------------------------------------------------------------
gitcmd(args...) = try readchomp(Cmd(`git -C $ROOT $args`)) catch; "unknown" end

function cc_flags()
    sh = joinpath(ROOT, "deps", "model", "compile.sh")
    m = isfile(sh) ? match(r"gcc -c ([^\n]*?) \w+\.c", read(sh, String)) : nothing
    return m === nothing ? "unknown" : m.captures[1]
end

function metadata(opts)
    cc_version = try first(split(readchomp(`cc --version`), '\n')) catch; "unknown" end
    return Dict{String,Any}(
        "timestamp" => string(now()),
        "tag" => opts["tag"],
        "git_commit" => gitcmd("rev-parse", "HEAD"),
        "git_branch" => gitcmd("rev-parse", "--abbrev-ref", "HEAD"),
        "git_dirty" => !isempty(gitcmd("status", "--porcelain", "--untracked-files=no")),
        "julia_version" => string(VERSION),
        "julia_threads" => Threads.nthreads(),
        "cpu" => Sys.cpu_info()[1].model,
        "cpu_threads" => Sys.CPU_THREADS,
        "os" => string(Sys.KERNEL, " ", Sys.MACHINE),
        "hostname" => gethostname(),
        "c_compiler" => cc_version,
        "c_flags" => cc_flags(),
        "quick" => opts["quick"],
        "filter" => something(opts["filter"], ""),
    )
end

# ------------------------------------------------------------------------------------------
# Results I/O
# ------------------------------------------------------------------------------------------
function save_run(results, meta)
    sha = meta["git_commit"][1:min(8, end)] * (meta["git_dirty"] ? "+dirty" : "")
    name = Dates.format(now(), "yyyymmdd-HHMMSS") * "_" * sha * (isempty(meta["tag"]) ? "" : "_" * meta["tag"])
    dir = joinpath(RESULTS_DIR, name)
    mkpath(dir)
    BenchmarkTools.save(joinpath(dir, "results.json"), results)
    open(io -> TOML.print(io, meta), joinpath(dir, "meta.toml"), "w")
    return dir
end

function resolve_run(spec, exclude)
    isdir(spec) && return abspath(spec)
    isdir(RESULTS_DIR) || error("no saved runs in $RESULTS_DIR")
    runs = sort(filter(d -> isdir(joinpath(RESULTS_DIR, d)) && joinpath(RESULTS_DIR, d) != exclude,
                       readdir(RESULTS_DIR)); rev=true)
    if spec != "latest"
        runs = filter(d -> TOML.parsefile(joinpath(RESULTS_DIR, d, "meta.toml"))["tag"] == spec, runs)
    end
    isempty(runs) && error("no saved run matches \"$spec\"")
    return joinpath(RESULTS_DIR, first(runs))
end

load_run(dir) = (BenchmarkTools.load(joinpath(dir, "results.json"))[1],
                 TOML.parsefile(joinpath(dir, "meta.toml")))

# ------------------------------------------------------------------------------------------
# Reporting
# ------------------------------------------------------------------------------------------
fmt_time(ns) = ns < 1e3 ? @sprintf("%.0f ns", ns) : ns < 1e6 ? @sprintf("%.1f µs", ns / 1e3) :
               ns < 1e9 ? @sprintf("%.2f ms", ns / 1e6) : @sprintf("%.3f s", ns / 1e9)
fmt_mem(b) = b < 1024 ? "$(b) B" : b < 1024^2 ? @sprintf("%.1f KiB", b / 1024) :
             b < 1024^3 ? @sprintf("%.1f MiB", b / 1024^2) : @sprintf("%.2f GiB", b / 1024^3)

# Find the trial at `path` in a (nested) BenchmarkGroup, or nothing
function lookup(g, path)
    for k in path
        (g isa BenchmarkGroup && haskey(g, k)) || return nothing
        g = g[k]
    end
    return g
end

function report(io, results, base; tol=0.05)
    rows = sort(collect(BenchmarkTools.leaves(results)); by=first)
    head = ["benchmark", "min", "median", "cv", "allocs", "memory", "ns/ch·smp"]
    base === nothing || append!(head, ["base median", "ratio", "verdict"])
    table = Vector{Vector{String}}()
    for (path, trial) in rows
        key = join(path, "/")
        med, mn = median(trial), minimum(trial)
        cv = length(trial.times) > 1 ? std(trial.times) / mean(trial.times) : NaN
        work = get(WORK, key, nothing)
        row = [key, fmt_time(time(mn)), fmt_time(time(med)), isnan(cv) ? "-" : @sprintf("%.1f%%", 100cv),
               string(allocs(med)), fmt_mem(memory(med)),
               work === nothing ? "-" : @sprintf("%.1f", time(med) / work)]
        if base !== nothing
            bt = lookup(base, path)
            if bt !== nothing
                b = median(bt)
                j = judge(med, b; time_tolerance=tol)
                v = BenchmarkTools.time(j)
                push!(row, fmt_time(time(b)), @sprintf("%.3f", time(med) / time(b)),
                      v == :improvement ? "faster" : v == :regression ? "SLOWER" : "~")
            else
                push!(row, "-", "-", "new")
            end
        end
        push!(table, row)
    end
    widths = [maximum(length, [head[k]; [r[k] for r in table]]) for k in eachindex(head)]
    line(cells) = join([k == 1 ? rpad(c, widths[k]) : lpad(c, widths[k]) for (k, c) in enumerate(cells)], "  ")
    println(io, line(head))
    println(io, "-"^sum(widths .+ 2))
    foreach(r -> println(io, line(r)), table)
end

# ------------------------------------------------------------------------------------------
# Main
# ------------------------------------------------------------------------------------------
"""
    main(args::Vector{String}=String[])

Run the benchmark suite with command-line style `args` (see the header of this file).
Returns the `BenchmarkGroup` of results, or `nothing` for `--list`.
"""
function main(args::AbstractVector{<:AbstractString}=String[])
    OPTS = parse_args(args)

    suite = SUITE
    if OPTS["filter"] !== nothing
        re = Regex(OPTS["filter"])
        suite = BenchmarkGroup()
        for (p, b) in BenchmarkTools.leaves(SUITE)
            occursin(re, join(p, "/")) || continue
            haskey(suite, p[1]) || (suite[p[1]] = BenchmarkGroup())
            sub = suite[p[1]]
            sub[p[2]] = b
        end
    end

    if OPTS["list"]
        foreach(((p, _),) -> println(join(p, "/")), sort(collect(BenchmarkTools.leaves(suite)); by=first))
        return nothing
    end

    n_bench = length(collect(BenchmarkTools.leaves(suite)))
    n_bench == 0 && error("no benchmarks match the filter")
    if OPTS["quick"]
        foreach(((_, b),) -> (b.params.seconds = 0.5), BenchmarkTools.leaves(suite))
    end

    meta = metadata(OPTS)
    println("Benchmarking $n_bench benchmarks | commit $(meta["git_commit"][1:min(8, end)])",
            meta["git_dirty"] ? " (dirty)" : "", " on $(meta["git_branch"]) | Julia $(meta["julia_version"]) ",
            "($(meta["julia_threads"]) threads) | $(meta["cpu"]) | C flags: $(meta["c_flags"])")
    OPTS["quick"] && println("(quick mode: numbers are noisy; do not use for decisions)")

    t0 = time()
    results = run(suite; verbose=false)
    println("Finished in $(round(time() - t0; digits=1)) s\n")

    base = nothing
    if OPTS["compare"] !== nothing
        base_dir = resolve_run(OPTS["compare"], "")
        base, base_meta = load_run(base_dir)
        println("Comparing against $(basename(base_dir)) (tag: \"$(base_meta["tag"])\"; ",
                "commit $(base_meta["git_commit"][1:min(8, end)])$(base_meta["git_dirty"] ? "+dirty" : ""))")
        if base_meta["cpu"] != meta["cpu"] || base_meta["hostname"] != meta["hostname"]
            println("WARNING: baseline was recorded on a different machine ($(base_meta["hostname"]), $(base_meta["cpu"]))")
        end
        println("Verdicts compare medians with a ±$(round(Int, 100 * OPTS["tol"]))% tolerance\n")
    end

    report(stdout, results, base; tol=OPTS["tol"])

    if OPTS["save"]
        dir = save_run(results, meta)
        println("\nSaved to $(relpath(dir, ROOT))")
    end
    return results
end

end # module BenchmarkRunner

if abspath(PROGRAM_FILE) == @__FILE__
    BenchmarkRunner.main(ARGS)
end
