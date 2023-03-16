using Helios
using AuditorySignalUtils
using AuditoryNerveFiber
using DrWatson
using CairoMakie
using DSP
update_theme!(fontsize=20)

# Synthesize pure-tone stimulus
x = scale_dbspl(cosine_ramp(pure_tone(1000.0, 0.0, 0.2, 100e3), 0.01, 100e3), 50.0)

# Simulate HSR/LSR/CN/IC responses
resp = sim_gfc2023_dict(
    vcat(x, zeros(10000)), 
    1000.0; 
    cn_tau_e=1.0e-3,
    cn_tau_i=1.0e-3,
    cn_delay=1.0e-3,
    cn_amp=1.0,
    cn_inh=0.5,
    ic_tau_e=1.0e-3,
    ic_tau_i=2.0e-3,
    ic_delay=2.0e-3,
    ic_amp=1.0,
    ic_inh=0.8,
    moc_cutoff=0.5,
    moc_beta=0.010,
    moc_weight_wdr=0.0,
    moc_weight_ic=0.0,
    dur_pad_left=0.05,
);
orig = sim_orig_dict(
    vcat(x, zeros(10000)),
    1000.0;
)

# Plot both responses
fig = Figure()
axs = map(1:4) do idx
    Axis(fig[idx, 1])
end
map(zip(axs, [resp["hsr"], resp["lsr"], resp["cn"], resp["ic"]])) do (ax, r)
    lines!(ax, r)
end
ylims!(axs[1], 0.0, 1000.0)
ylims!(axs[2], 0.0, 200.0)
ylims!(axs[3], 0.0, 220.0)
ylims!(axs[4], 0.0, 62.0)
#lines!(axs[3], 1.2 .* sim_sfie_nc2004(resp["hsr"]; τ_e=1e-3, τ_i=1e-3, d_i=1e-3, S=0.5, A=1.0); color=:orange)
#lines!(axs[3], 1.2 .* sim_sfie_nc2004(orig["hsr"]; τ_e=1e-3, τ_i=1e-3, d_i=1e-3, S=0.5, A=1.0); color=:purple)
fig

# Simpler version
fig = Figure()
ax = Axis(fig[1, 1])
x = vcat(1.0, zeros(9999))
b, a = AuditoryMidbrain.get_α_normalized(1e-3, 100e3, 1.0)
E = filt(b, a, x)
I = filt(b, a, shiftsignal(x, Int(floor(d_i*fs))))
lines!(ax, E; color=:red)
lines!(ax, I; color=:blue)
ylims!(ax, 0.0, 400)
xlims!(ax, 0.0, 1000.0)
fig