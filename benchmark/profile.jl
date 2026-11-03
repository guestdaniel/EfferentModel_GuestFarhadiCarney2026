# Statistical profile of a representative workload, with C-level attribution.
#
#   julia --project=benchmark benchmark/profile.jl [options]   # standalone, from the shell
#
#   include("benchmark/profile.jl")                            # from a persistent session
#   ProfileRunner.main(["--chan", "32", "--dur", "0.5"])       # (e.g. the Julia MCP server)
#
# Options:
#   --chan N        Number of CF channels (default 32)
#   --dur S         Stimulus duration in seconds (default 0.5)
#   --binaural      Profile the binaural model
#   --kw 'EXPR'     Extra keyword arguments for sim_gfc2023, e.g. --kw 'moc_width=0.0'
#   --reps N        Number of repetitions to profile (default 3)
#   --top N         Rows to show in each table (default 25)
#
# Reports (1) self time ("leaf" samples) by function and by category (model C code, libm,
# GC, memory, Julia wrapper, ...) and (2) inclusive time of exported functions in the model
# library. Only awake samples from the thread running the workload are counted, so idle
# threads and parallel GC worker threads are excluded.
#
# Note: helper functions in model.c (C1ChirpFilt, WbGammaTone, ...) are usually inlined into
# model_monaural/model_binaural at -O3, so their self time is reported under those names.
# libm calls (pow, exp, atan, ...) are never inlined and are attributed exactly.

module ProfileRunner

using Profile
using Printf
using Libdl
using caplab_binaural_efferent_model
const BEM = caplab_binaural_efferent_model

const FS = 100e3

const MODEL_SYMS = Set(["model_monaural", "model_binaural", "middle_ear", "C1ChirpFilt",
    "C2ChirpFilt", "WbGammaTone", "OhcLowPass", "IhcLowPass", "NLafterohc", "NLogarithm",
    "Boltzman", "gain_groupdelay", "Get_taubm", "Get_tauwb", "apply_ws1988_adaptation",
    "apply_powerlaw_adaptation", "apply_powerlaw_adaptation_iir", "filter_alpha",
    "get_alpha_norm", "moc_nonlinearity", "normal_pdf", "delay_signal", "filter_lowpass_iir",
    "hw_rectify"])

function parse_opts(args)
    o = Dict{String,Any}("chan" => 32, "dur" => 0.5, "binaural" => false, "kw" => "", "reps" => 3, "top" => 25)
    i = 1
    while i <= length(args)
        a = args[i]
        if a == "--binaural"
            o["binaural"] = true; i += 1
        elseif a in ("--chan", "--reps", "--top")
            o[a[3:end]] = parse(Int, args[i+1]); i += 2
        elseif a == "--dur"
            o["dur"] = parse(Float64, args[i+1]); i += 2
        elseif a == "--kw"
            o["kw"] = args[i+1]; i += 2
        else
            error("unknown argument: $a")
        end
    end
    return o
end

# Named function so that its frames can be found in the profile data
function workload(x, cf, kw, reps)
    for _ in 1:reps
        sim_gfc2023(x, cf; kw...)
    end
end

libname(fr) = basename(String(fr.file))

function category(fr)
    f, file = String(fr.func), String(fr.file)
    fr.from_c || return occursin("caplab_binaural_efferent_model", file) ? "Julia: package wrappers" : "Julia: Base/other"
    occursin("libsystem_m", file) && return "libm (transcendentals)"
    (occursin("gc", file) || startswith(f, "jl_gc") || startswith(f, "ijl_gc") || startswith(f, "gc_")) && return "GC"
    (occursin("libsystem_platform", file) || f in ("memset", "memmove", "memcpy", "bzero") ||
        startswith(f, "_platform_mem")) && return "memset/memcpy"
    occursin("libsystem_malloc", file) && return "malloc/free"
    (occursin("dyld", file) || f in ("dlopen", "dlclose", "dlsym")) && return "dynamic loading"
    occursin("libgfc2026", file) && return "model C code"
    return "other C"
end

function print_table(title, d, n, top)
    println("\n", title)
    for (k, v) in first(sort(collect(d); by=last, rev=true), top)
        println("  ", @sprintf("%5.1f%%", 100v / n), "  ", k)
    end
end

"""
    main(args::Vector{String}=String[])

Profile `sim_gfc2023` with command-line style `args` (see the header of this file) and
print the attribution tables.
"""
function main(args::AbstractVector{<:AbstractString}=String[])
    O = parse_opts(args)
    x1 = 0.02 .* randn(round(Int, O["dur"] * FS))
    x = O["binaural"] ? [x1, 0.02 .* randn(length(x1))] : x1
    cf = O["chan"] == 1 ? [1000.0] : collect(exp.(range(log(125.0), log(16e3), O["chan"])))
    kw = isempty(O["kw"]) ? NamedTuple() : Core.eval(@__MODULE__, Meta.parse("(; $(O["kw"]))"))

    # Keep the model library loaded for the whole run so that its addresses can be
    # symbolized (model! dlopens/dlcloses it on every call, which would otherwise unload it)
    lib = Libdl.dlopen(BEM.libgfc2026)
    try
        sim_gfc2023(x, cf; kw...)                       # compile first
        Profile.init(n=10^7, delay=0.0002)
        Profile.clear()
        t = @elapsed @profile workload(x, cf, kw, O["reps"])
        data = Profile.fetch()   # with metadata (thread id, sleep state) per sample
        lidict = Profile.getdict(data)

        # Split raw data into per-sample blocks: [ips..., threadid, taskid, cycles, sleepstate, 0, 0]
        blocks = Tuple{Vector{UInt64},UInt64,UInt64}[]   # (backtrace, threadid, sleepstate)
        start = 1
        for i in eachindex(data)
            if Profile.is_block_end(data, i)
                push!(blocks, (data[start:i-6], data[i-Profile.META_OFFSET_THREADID],
                               data[i-Profile.META_OFFSET_SLEEPSTATE]))
                start = i + 1
            end
        end
        frames(bt) = Iterators.flatten(lidict[ip] for ip in bt)

        # Identify the thread running the workload, then keep all of its awake samples.
        # (Filtering on "stack contains the workload" instead would drop samples in libm
        # leaf functions, whose stacks often fail to unwind past the leaf.)
        wl_threads = [tid for (bt, tid, _) in blocks if any(fr -> fr.func === :workload, frames(bt))]
        isempty(wl_threads) && error("no workload samples collected")
        wl_tid = argmax(t -> count(==(t), wl_threads), unique(wl_threads))
        samples = [bt for (bt, tid, sleep) in blocks if tid == wl_tid && sleep == 1 && !isempty(bt)]
        n = length(samples)
        n == 0 && error("no workload samples collected")

        self_fn = Dict{String,Int}()
        self_cat = Dict{String,Int}()
        incl = Dict{String,Int}()
        for bt in samples
            fr = first(lidict[first(bt)])
            name = fr.from_c ? "$(fr.func) [$(libname(fr))]" : "$(fr.func) [$(libname(fr)):$(fr.line)]"
            self_fn[name] = get(self_fn, name, 0) + 1
            c = category(fr)
            self_cat[c] = get(self_cat, c, 0) + 1
            for f in unique(String(fr2.func) for fr2 in frames(bt) if String(fr2.func) in MODEL_SYMS)
                incl[f] = get(incl, f, 0) + 1
            end
        end

        println("Workload: $(O["binaural"] ? "binaural" : "monaural"), $(O["chan"]) channels, ",
                "$(O["dur"]) s stimulus, kwargs: $(isempty(O["kw"]) ? "(defaults)" : O["kw"]), ",
                "$(O["reps"]) reps in $(round(t; digits=2)) s; $n samples")
        print_table("Self time by category", self_cat, n, O["top"])
        print_table("Self time by function (leaf frame)", self_fn, n, O["top"])
        print_table("Inclusive time, model library functions (not inlined)", incl, n, O["top"])
    finally
        Libdl.dlclose(lib)
    end
    return nothing
end

end # module ProfileRunner

if abspath(PROGRAM_FILE) == @__FILE__
    ProfileRunner.main(ARGS)
end
