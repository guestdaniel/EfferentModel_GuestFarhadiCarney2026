# This script provides the standard profiling and benchmarking scripts for the efferent
# model code. The broader goal of this script is to identify major inefficiencies or low-
# hanging fruit for optimization within the model C code and/or Julia/MATLAB wrappers.
# We begin at the lowest possible level, comparing simple model function calls to a matched
# benchmark written entirely in C with minimal boilerplate. Then we examine the behavior of
# the model as we make many repeated function calls (high n_rep) or simulate many channels
# (high n_chan) to look for systematic issues with scaling.
using Helios
using Profile
using ProfileView
using BenchmarkTools

# Define function that represents a standard model function call.
# We use a 60 dB sinusoid as input, default 5 seconds with one channel. This is a nice
# analogue to the test script defined in test.c. Note that we set moc_weight > 0 and
# moc_offset > 0 to simulate normal conditions in terms of how often we need to compute moc
# nonlinearities.
function test(; N=500_000, n_chan=20)
    stim = 10.0 ^ (60.0 / 20.0) * 20e-6 .* sin.(2.0 .* π .* ((0:N) ./ 100e3))
    sim_gfc2023_dict(stim, fill(1e3, n_chan); moc_weight=fill(10.0, n_chan), moc_offset=1.0, dur_pad_left=0.0, dur_pad_right=0.0)
end

# Observation:
#   test.c program @ 5 seconds one channel with -o3 => ~1.15 seconds
#   matching @benchmark below =>                       ~1.12 seconds???
#
#   Ergo, Julia code for single-channel runs has minimal overhead, except when GC has
#   to occur. This may scale poorly in multi-channel run cases. Profiling code below
#   additionally supports this interpretation (relatively little of total execution time
#   is spent anywhere but libgfc2023.so!)

# Benchmark 5-s one-channel run to compare to available direct C test script
@benchmark test(; N=500_000, n_chan=1) samples=50 

# Profile same model function call
Profile.clear()
@profile test(; N=500_000, n_chan=1)
ProfileView.view(; C=true)
# Profile.print(; format=:flat, sortedby=:count, C=true)

# Compare above case to a multi-channel case; does the ratio of overhead to compute stay constant?
Profile.clear()
@profile test(; N=100_000, n_chan=5)
ProfileView.view(; C=true)

# Observation:
#   Model initially scales well with increasing number of channels but eventually
#   slows down relative to predictions from n_chan=1 baseline. I presume this reflects the
#   accumulated effects of garbage collection in some way? Another possibility is accumulated
#   penalties associated with cache misses.

# Plot growth with number of channels
n_chans = 2 .^ (0:1:7)
n_rep = 10 
vals = map(n_chans) do n_chan
    map(1:n_rep) do _
        @elapsed test(; N=10_000, n_chan=n_chan)
    end
end
fig = Figure()
ax = Axis(fig[1, 1]; xlabel="Num chans", ylabel="Runtime (s)")
lines!(ax, n_chans, mean.(vals); color=:black)
errorbars!(ax, n_chans, mean.(vals), std.(vals); color=:black, whiskerwidth=5.0)
lines!(ax, n_chans, n_chans .* mean(vals[1]); color=:gray, linestyle=:dash)
ylims!(ax, 0.0, maximum(mean.(vals)) * 1.2)
xlims!(ax, 0.0, maximum(n_chans)+1)
fig

# Observation:
#   

# Benchmark short multichannel vs long single channel
@benchmark test(; N=100_000, n_chan=5) samples=20
@benchmark test(; N=500_000, n_chan=1) samples=20

