using CairoMakie
using Helios_plots
using Helios
using CarneyLabUtils2
using AuditorySignalUtils
using DSP

### ========================================================================================
# Look at SAM tone, 8 kHz, 40 Hz mod rate, moderate sound level, etc.
### ========================================================================================
cf_center = 8000.0
cf_range = 2/3
n_cf = 29
cf = LogRange(cf_center * 2 ^ (-cf_range), cf_center * 2 ^ cf_range, n_cf)
x = sam(8000.0, 80.0)
#x = sam(8000.0, 40.0)
resp = sim_gfc2023_dict(
    x, 
    cf; 
    dur_pad_right=0.02, 
    clip_right=false, 
    moc_cutoff=0.2,
    moc_weight_wdr=3.0,
    moc_weight_ic=2.0,
    # moc_weight_wdr=0.0,
    # moc_weight_ic=0.0,
    ic_amp=2.0,
)

# Plot time course at CF
plot_timecourse!(resp, Int(ceil(n_cf/2)), ["ihc", "hsr", "lsr", "ic"], Figure(; resolution=(1200, 1100)))

# Plot stacks of timecourses across CF
plot_timecourse_stack!(cf, resp, ["hsr", "lsr", "ic"], 3:3:29, Figure(; resolution=(1600, 2000)))

# Plot neurograms
neurograms = map(["hsr", "lsr", "ic", "gain"]) do output
    fig, ax = plot_neurogram(cf, resp[output])
    return fig
end
comb = displayimg(hcat(getimg.(neurograms)...))

### ========================================================================================
# Look at two-formant stimulus
### ========================================================================================
n_cf = 41
cf = LogRange(800.0, 8000.0, n_cf)
fs = 100e3
f0 = 100.0

function two_pole_resonator(x, ω, r)
    a = [1.0, -2*r*cos(ω), r^2]
    b = [1.0]
    filt(b, a, x)
end

function two_pole_resonator(x, f, bw, fs)
    ω = 2π * (f/fs)
    r = exp(-π * bw/fs)
    two_pole_resonator(x, ω, r)
end


# Synthesize two-formant complex tone
x = sum(map(f -> cosine_ramp(scale_dbspl(pure_tone(f, 0.0, 0.3, fs), 50.0), 0.01, fs), f0:f0:min(10e3, fs/2)))
x = two_pole_resonator(x, 1500.0, 400.0, fs) .+ 2.0 .* two_pole_resonator(x, 3000.0, 400.0, fs)
x = scale_dbspl(x, 70.0)

@time resp = sim_gfc2023_dict(
    x, 
    cf; 
    dur_pad_right=0.02, 
    clip_right=false, 
    moc_cutoff=0.2,
    moc_weight_wdr=3.0,
    moc_weight_ic=1.0,
    # moc_weight_wdr=0.0,
    # moc_weight_ic=0.0,
    ic_amp=2.0,
)

# Plot time course at CF
plot_timecourse!(resp, Int(ceil(n_cf/2)), ["ihc", "hsr", "lsr", "ic"], Figure(; resolution=(1200, 1100)))

# Plot stacks of timecourses across CF
plot_timecourse_stack!(cf, resp, ["hsr", "lsr", "ic"], 3:3:n_cf, Figure(; resolution=(1600, 2000)))

# Plot neurograms
neurograms = map(["hsr", "lsr", "ic", "gain"]) do output
    fig=Figure(; resolution=(1800, 800))
    fig, ax = plot_neurogram(x, cf, resp[output], fig)
    xlims!(ax, 0.05, 0.15)
    hlines!(ax, [1500.0, 3000.0]; color=:white, style=:dash)
    return fig
end
comb = displayimg(hcat(getimg.(neurograms)...))

# Plot time courses at channels closest to peaks
chans = map(x -> argmin(abs.(cf .- x)), [1500.0, 3000.0])
plot_timecourse_stack!(cf, resp, ["hsr", "lsr", "ic", "mocwdr", "mocic"], chans, Figure(; resolution=(2000, 700)); title=true)
