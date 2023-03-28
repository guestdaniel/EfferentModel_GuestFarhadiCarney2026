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

