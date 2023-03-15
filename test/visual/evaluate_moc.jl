using Helios
using AuditorySignalUtils
using AuditoryNerveFiber
using DrWatson
using CairoMakie
using DSP
pt(f=1000.0, l=50.0, dur=0.2, fs=100e3) = scale_dbspl(pure_tone(f, 0.0, dur, fs), l)
update_theme!(fontsize=20)

### ========================================================================================
### Tone responses
### ========================================================================================
# Simulate MOC response to 1-kHz pure tone
resp = sim_gfc2023_dict(pt(1000.0, 50.0, 0.2), 1000.0; moc_cutoff=10.0);

# Plot response
fig = Figure()
ax = Axis(fig[1, 1])
t = 0.0:(1/100e3):(0.2-1/100e3)
lines!(ax, t, resp["lsr"]; color=:black)
lines!(ax, t, resp["moc"]; color=:cyan)
#xlims!(ax, 0.01, 0.015)
fig