using CairoMakie
using Helios_plots
using Helios
using CarneyLabUtils2

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
    ic_amp=2.0,
)
plot_timecourse!(resp, ["stim", "ihc", "hsr", "lsr", "ic"], Figure(; resolution=(1200, 1100)))

# Look at SAM tones
x = sam(8000.0, 80.0)
resp = sim_gfc2023_dict(
    x, 
    8000.0; 
    dur_pad_right=0.02, 
    clip_right=false, 
    moc_cutoff=0.2,
    moc_weight_wdr=2.0,
    moc_weight_ic=2.0,
    ic_amp=2.0,
)
plot_timecourse!(resp, ["stim", "ihc", "hsr", "lsr", "ic"], Figure(; resolution=(1200, 1100)))

# Look at MTF with and without efferent feedback
m1 = AuditorySubcortexGFCv1(; moc_weight_wdr=0.0, moc_weight_ic=0.0)
m2 = AuditorySubcortexGFCv1(; moc_weight_wdr=2.0, moc_weight_ic=2.0)
m3 = AuditorySubcortexGFCv1(; moc_weight_wdr=2.0, moc_weight_ic=4.0)
analyses = map([m1, m2, m3]) do m
    analysis = AnalysisNoiseMTFModel(m)
    load_or_run!(analysis)
    return analysis
end

viz(analyses[1].mtf)
viz(analyses[2].mtf)