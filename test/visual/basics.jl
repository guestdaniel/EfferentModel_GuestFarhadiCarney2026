using Helios
using AuditorySignalUtils
using AuditoryNerveFiber
using DrWatson
using CairoMakie
using DSP
using Statistics

pt(f=1000.0, l=50.0, dur=0.2, fs=100e3) = scale_dbspl(pure_tone(f, 0.0, dur, fs), l)

function iso_level_curve(; level=50.0, freqs=LogRange(200.0, 20e3, 40))
    # Map through frequencies, synthesize tone and compute mean response
    rates = map(freqs) do freq
        stim = pt(freq, level)
        mean(sim_anrate_gfc2023(stim, 1000.0))
    end 

    # Plot
    fig = Figure()
    ax = Axis(fig[1, 1]; xscale=log10)
    lines!(ax, freqs, rates)
    ax.xlabel = "Frequency (Hz)"
    ax.ylabel = "Rate (sp/s)"
    fig
end
fig = iso_level_curve()
save(projectdir("outputs", "test", "iso_level_curve.png"), fig)

function excitation_pattern(; freq=1000.0, level=50.0, cfs=LogRange(200.0, 20e3, 40))
    # Synthesize tone
    x = pt(freq, level)

    # Simulate response at all CFs simultaneously
    rates = mean.(sim_anrate_gfc2023(x, cfs))

    # Plot
    fig = Figure()
    ax = Axis(fig[1, 1]; xscale=log10)
    lines!(ax, cfs, rates)
    ax.xlabel = "CF (Hz)"
    ax.ylabel = "Rate (sp/s)"
    fig
end
fig = excitation_pattern()
save(projectdir("outputs", "test", "excitation_pattern.png"), fig)
