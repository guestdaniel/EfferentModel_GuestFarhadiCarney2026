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
resp = sim_gfc2023_dict(
    pt(1000.0, 50.0, 0.2), 
#    zeros(20_000),
    1000.0; 
    moc_cutoff=1.0,
    moc_beta=0.01,
    moc_offset=0.0,
    moc_maxrate=50.0,
);

# Plot response
fig = Figure()
ax = Axis(fig[1, 1])
t = 0.0:(1/100e3):(0.2-1/100e3)
lines!(ax, t, resp["lsr"]; color=:black)
lines!(ax, t, resp["wdr"]; color=:cyan)
lines!(ax, t, resp["moc"]; color=:pink)
fig

### ========================================================================================
### SAM responses
### ========================================================================================
# Simulate MOC response to 1-kHz SAM tone
resp = sim_gfc2023_dict(
    sam(1000.0, 20.0, -0.0, 50.0, 0.3), 
#    zeros(20_000),
    1000.0; 
    moc_cutoff=1.0,
    moc_beta=0.01,
    moc_offset=0.0,
    moc_maxrate=50.0,
);

# Plot response
fig = Figure()
ax = Axis(fig[1, 1])
t = 0.0:(1/100e3):(0.3-1/100e3)
lines!(ax, t, resp["lsr"]; color=:black)
lines!(ax, t, resp["wdr"]; color=:cyan)
lines!(ax, t, resp["moc"]; color=:pink)
fig