using Helios
using AuditorySignalUtils
using AuditoryNerveFiber
using DrWatson
using CairoMakie
using DSP
update_theme!(fontsize=20)

### ========================================================================================
### Plotting functions
### ========================================================================================
function plot_moc(resp)
    # Derive duration and time axis
    dur = length(resp["ihc"])/100e3
    t = 0.0:(1/100e3):nextfloat(dur-1/100e3)

    # Create plot
    fig = Figure()
    ax_wdr = Axis(fig[1, 1])
    ax_wdr_r = Axis(fig[1, 1]; yaxisposition=:right)
    ax_ic = Axis(fig[2, 1])
    ax_ic_r = Axis(fig[2, 1]; yaxisposition=:right)
    hidespines!.([ax_wdr_r, ax_ic_r])
    hidexdecorations!.([ax_wdr, ax_ic], ticks=false, ticklabels=false, label=false)
    hidexdecorations!(ax_wdr, ticks=false, label=false)
    hidexdecorations!.([ax_wdr_r, ax_ic_r])
    hideydecorations!.([ax_wdr, ax_wdr_r, ax_ic, ax_ic_r], ticks=false, ticklabels=false)
    ylims!.([ax_wdr_r, ax_ic_r], 0.0, 1.02)
    ylims!.(ax_wdr, 0.0, 1.2 * maximum(resp["lsr"]))
    ylims!.(ax_ic, 0.0, 1.2 * maximum(resp["ic"]))
    xlims!.([ax_wdr, ax_wdr_r, ax_ic, ax_ic_r], -0.01, maximum(t)*1.25)

    # Plot inputs
    lines!(ax_wdr, t, resp["lsr"]; color=:black)
    lines!(ax_ic, t, resp["ic"]; color=:black)

    # Plot lowpassed inputs
    lines!(ax_wdr, t, resp["wdr"]; color=:cyan)
    lines!(ax_ic, t, resp["icin"]; color=:cyan)

    # Plot gain
    lines!(ax_wdr_r, t, resp["moc"]; color=:pink, linewidth=3.0)
    lines!(ax_ic_r, t, resp["moc"]; color=:pink, linewidth=3.0)

    # Plot labels
    text!(ax_wdr_r, maximum(t)+0.01, resp["moc"][end]; text="MOC output", color=:pink, alignment=(:left, :center))
    text!(ax_ic_r, maximum(t)+0.01, resp["moc"][end]; text="MOC output", color=:pink, alignment=(:left, :center))
    text!(ax_wdr, maximum(t)+0.01, resp["wdr"][end]; text="MOC input", color=:cyan, alignment=(:left, :center))
    text!(ax_ic, maximum(t)+0.01, resp["icin"][end]; text="MOC input", color=:cyan, alignment=(:left, :center))

    # Adjust labels and finer details
    ax_ic.xlabel = "Time re: simulation start (s)"
    ax_wdr.ylabel = "Firing rate (sp/s)"
    ax_ic.ylabel = "Firing rate (sp/s)"
    ax_wdr_r.ylabel = "Cochlear gain (0, 1]"
    ax_ic_r.ylabel = "Cochlear gain (0, 1]"
    fig
end

### ========================================================================================
### Tone responses
### ========================================================================================
resp = sim_gfc2023_dict(
    pt(1000.0, 50.0, 0.2), 
    1000.0; 
    ic_tau_e=1.0e-3,
    ic_tau_i=2.0e-3,
    ic_delay=2.0e-3,
    ic_amp=1.0,
    ic_inh=0.9,
    moc_cutoff=0.5,
    moc_beta=0.010,
    moc_weight_wdr=3.0,
    moc_weight_ic=2.0,
    dur_pad_left=0.2,
    dur_pad_right=0.05,
    clip_left=false,
    clip_right=false,
);

# Plot response
plot_moc(resp)

### ========================================================================================
### SAM responses
### ========================================================================================
# Simulate MOC response to 1-kHz SAM tone
resp = sim_gfc2023_dict(
    vcat(sam(8000.0, 70.0, -6.0, 50.0, 0.3), zeros(10000)), 
    8000.0; 
    ic_tau_e=1.0e-3,
    ic_tau_i=2.0e-3,
    ic_delay=1.0e-3,
    ic_amp=0.2,
    ic_inh=0.8,
    moc_cutoff=0.5,
    moc_beta=0.010,
    moc_offset=10.0,
    moc_weight_wdr=5.0,
    moc_weight_ic=20.0,
    dur_pad_left=0.2,
    dur_pad_right=0.05,
    clip_left=false,
    clip_right=false,
);

# Plot response
plot_moc(resp)
