# profile_single_vs_repeated_call.jl
#
# This script compares runtime for the Helios model for a single long-duration simulation
# versus multiple short-duration simulations with the same overall time.
using Helios
using Profile
using ProfileView
using BenchmarkTools
using CairoMakie
using AuditoryNerveFiber

# Profile overall model call
function test(n=100_000, ratio=1, mode="old")
    # Assert that n can be divided evenly into ratio parts
    @assert isinteger(n/ratio) 

    # Run the model ratio times
    for ii = 1:ratio
        # Run the model for n/ratio length stimulus
        n_reduced = Int(n/ratio)

        # Run model
        if mode == "old"
            sim_ihc_zbc2014(zeros(n_reduced), 1e3)
            sim_anrate_zbc2014(zeros(n_reduced), 1e3)
            sim_anrate_zbc2014(zeros(n_reduced), 1e3)
        else
            sim_gfc2023(zeros(n_reduced), 1e3; dur_pad_right=0.0, dur_pad_left=0.0)
        end
    end
end

# Plot general result
fig = Figure()
ax = Axis(fig[1, 1]; yscale=log10, xscale=log10)
ratios = [1, 2, 4, 10, 20, 40, 100]
baselines = [100_000]
map(["old", "new"]) do m
    runtimes = map(ratios) do r
        return @elapsed test(100_000, r, m)
    end
    lines!(ax, ratios, runtimes; label=m)
end
axislegend()
ax.xticks = ratios
ylims!(ax, 0.01, 100.0)
fig

# Plot general result
fig = Figure()
ax = Axis(fig[1, 1]; yscale=log10, xscale=log10)
ratios = [1, 2, 4, 10, 20, 40, 100]
baselines = [100_000]
map(["old", "new"]) do m
    runtimes = map(ratios) do r
        return @elapsed test(100_000, r, m)
    end
    lines!(ax, ratios, runtimes ./ 100_000; label=m)
end
axislegend()
ax.xticks = ratios
ax.ylabel = "s/samp"
#ylims!(ax, 0.01, 100.0)
fig


# Benchmark 
@benchmark test(100_000, 1) samples=10
@benchmark test(100_000, 10) samples=10

# Profile
Profile.clear()
@profile test(1_000_000, 1)
#ProfileView.view(; C=true)
Profile.print(; format=:flat, sortedby=:count)

Profile.clear()
@profile test(1_000_000, 20)
#ProfileView.view(; C=true)
Profile.print(; format=:flat, sortedby=:count)
