using CairoMakie
using Helios_plots
using Helios
using CarneyLabUtils2
using AuditorySignalUtils
using DSP

### ========================================================================================
# Look at two-formant stimulus
### ========================================================================================
# Set parameters for stimulus and model
n_cf = 41
cf = LogRange(800.0, 8000.0, n_cf)
fs = 100e3
f0 = 100.0

# Synthesize two-formant complex tone
s = TwoFormantTone()
m = AuditorySubcortexGFCv1(;
    cf=cf,
    moc_cutoff=0.2,
    moc_weight_wdr=3.0,
    moc_weight_ic=1.0,
    ic_a=2.0,
    fractional=false,
    pathhead="\\home\\daniel\\cl_sim\\efferent\\explore_3-27"
)
if isfile(filename([m, s]))
    temp = load(filename([m, s]))
    x, resp = temp["x"], temp["resp"]
else
    x = synthesize!(s)
    resp = m(x)
    save(filename([m, s]), Dict("x" => x, "resp" => resp))
end

# Plot time course at CF
plot_timecourse!(resp, Int(ceil(n_cf/2)), ["ihc", "hsr", "lsr", "ic"], Figure(; resolution=(1200, 1100)))

# Plot stacks of timecourses across CF
plot_timecourse_stack!(cf, resp, ["hsr", "lsr", "ic"], 3:3:n_cf, Figure(; resolution=(1600, 2000)))

# Plot time courses at channels closest to peaks
chans = map(x -> argmin(abs.(cf .- x)), [1500.0, 3000.0])
plot_timecourse_stack!(cf, resp, ["hsr", "lsr", "ic", "mocwdr", "mocic"], chans, Figure(; resolution=(2000, 700)); title=true)

# Plot neurograms
neurograms = map(["ihc", "hsr", "lsr", "ic", "gain"]) do output
    fig=Figure(; resolution=(1800, 800))
    fig, ax = plot_neurogram(x, cf, resp[output], fig; clims=getylim(output), ylims=getylim_avg(output), cmap=getcmap(output))
    return fig
end
fig = displayimg(hcat(getimg.(neurograms)...))
save("\\home\\daniel\\cl_fig\\efferent\\3-27\\two_formant_neurograms.png", fig)

# Plot time courses at channels closest to peaks
plot_timecourse_stack!(cf, resp, ["hsr"], [12, 14], Figure(; resolution=(2000, 700)); title=true)
