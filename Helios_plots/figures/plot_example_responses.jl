using CairoMakie
using Helios_plots
using Helios

# Look at SAM tones
x = sam(8000.0, 40.0)
resp = sim_gfc2023_dict(
    x, 
    8000.0; 
    dur_pad_right=0.02, 
    clip_right=false, 
    moc_cutoff=0.2,
    moc_weight_wdr=2.0,
    moc_weight_ic=2.0,
)
plot_timecourse!(resp, ["stim", "ihc", "hsr", "lsr", "ic"], Figure(; resolution=(1200, 1100)))