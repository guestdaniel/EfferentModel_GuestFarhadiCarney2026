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
    t = 0.0:(1/100e3):(dur-1/100e3)

    # Create plot
    fig = Figure()
    ax_wdr = Axis(fig[1, 1])
    ax_wdr_r = Axis(fig[1, 1]; yaxisposition=:right)
    ax_ic = Axis(fig[2, 1])
    ax_ic_r = Axis(fig[2, 1]; yaxisposition=:right)
    hidespines!.([ax_wdr_r, ax_ic_r])
    hidexdecorations!.([ax_wdr, ax_ic], ticks=false, ticklabels=false)
    hidexdecorations!.([ax_wdr_r, ax_ic_r])
    hideydecorations!.([ax_wdr, ax_wdr_r, ax_ic, ax_ic_r], ticks=false, ticklabels=false)
#    linkaxes!(ax_wdr, ax_wdr_r)
#    linkaxes!(ax_ic, ax_ic_r)
    ylims!.([ax_wdr_r, ax_ic_r], 0.0, 1.0)
    ylims!.(ax_wdr, 0.0, 1.2 * maximum(resp["lsr"]))
    ylims!.(ax_ic, 0.0, 1.2 * maximum(resp["ic"]))

    # Plot inputs
    lines!(ax_wdr, t, resp["lsr"]; color=:black)
    lines!(ax_ic, t, resp["ic"]; color=:black)

    # Plot lowpassed inputs
    lines!(ax_wdr, t, resp["wdr"]; color=:cyan)
    lines!(ax_ic, t, resp["icin"]; color=:cyan)

    # Plot gain
    lines!(ax_wdr_r, t, resp["moc"]; color=:pink, linewidth=3.0)
    lines!(ax_ic_r, t, resp["moc"]; color=:pink, linewidth=3.0)
    fig
end

### ========================================================================================
### Tone responses
### ========================================================================================
resp = sim_gfc2023_dict(
    pt(1000.0, 50.0, 0.2), 
    1000.0; 
    moc_cutoff=1.0,
    moc_beta=0.01,
    moc_weight_wdr=1.0,
    moc_weight_ic=10.0,
);

# Plot response
plot_moc(resp)

### ========================================================================================
### SAM responses
### ========================================================================================
# Simulate MOC response to 1-kHz SAM tone
resp = sim_gfc2023_dict(
    sam(1000.0, 80.0, -0.0, 50.0, 0.3), 
    1000.0; 
    moc_cutoff=1.0,
    moc_beta=0.010,
    moc_weight_wdr=1.0,
    moc_weight_ic=3.0,
);

# Plot response
plot_moc(resp)
