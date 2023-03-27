using CairoMakie
using Helios
using CarneyLabUtils2
using AuditorySignalUtils
using AuditoryNerveFiber

# Simulate response to sam tones
x = sam(2000.0, 40.0)
cf = LogRange(2000.0 * 2^-0.5, 2000.0 * 2^0.5, 9)
resp = sim_gfc2023_dict(x, cf)
fig, ax = plot_neurogram(cf, hcat(resp["hsr"]...))
fig

x = sam(2000.0, 40.0)
cf = LogRange(2000.0 * 2^-0.5, 2000.0 * 2^0.5, 9)
orig = map(_cf -> sim_anrate_zbc2014(sim_ihc_zbc2014(x, _cf), _cf; fractional=false), cf)
fig, ax = plot_neurogram(cf, hcat(orig...))
fig
