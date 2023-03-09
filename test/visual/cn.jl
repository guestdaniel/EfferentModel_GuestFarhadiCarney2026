using Test
using AuditorySignalUtils
using AuditoryNerveFiber
using AuditoryMidbrain
using Statistics
using DSP
using Helios
using UnicodePlots

# Compare synapse outputs side by side 
pt(f=1000.0, l=50.0, dur=0.2, fs=100e3) = scale_dbspl(pure_tone(f, 0.0, dur, fs), l)
syn = sim_gfc2023_dict(pt(), 1000.0)["syn"];
cn_old = AuditoryMidbrain.sim_sfie_nc2004(syn)
cn_new = sim_gfc2023_dict(pt(), 1000.0)["cn"];
lineplot(cn_old[1:1000])
lineplot(cn_new[1:1000])
lineplot(cn_old .- cn_new)
