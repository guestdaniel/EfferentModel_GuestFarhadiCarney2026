using Helios
using AuditorySignalUtils
using AuditoryNerveFiber
using DrWatson
using CairoMakie

pt(f=1000.0, l=50.0, dur=0.2, fs=100e3) = scale_dbspl(pure_tone(f, 0.0, dur, fs), l)

function plot_model!()
    # Stimulus
    x = pt(1000.0, 50.0)

    # Responses
    resps_orig = sim_orig_dict(x, 1000.0)
    resps_new = sim_gfc2023_dict(x, 1000.0)
    stages = ["control", "c1", "c2", "ihc", "expon", "powerlaw"]

    # Plot output
    fig = Figure(; resolution=(800, 200*length(stages)))
    axs = [Axis(fig[i, 1]) for i in eachindex(stages)]
    for (ax, stage) in zip(axs, stages)
        lines!(ax, resps_orig[stage]; color=:gray)
        lines!(ax, resps_new[stage]; color=:red, linestyle=:dash, linewidth=3.0)
        ax.ylabel = stage
    end

    # Adjust labels
    xlims!.(axs, 5000, 6000)
    axs[end].xlabel = "Time (s)"
    hidexdecorations!.(axs[1:(end-1)], ticks=false, grid=false)
    fig
    save(projectdir("outputs", "test", "responses_up_to_powerlaw.png"), fig)
end
plot_model!()

function plot_model_2!()
    # Stimulus
    x = pt(1000.0, 50.0)

    # Responses
    resps_orig = sim_orig_dict(x, 1000.0)
    resps_new = sim_gfc2023_dict(x, 1000.0)
    stages = ["sout1", "sout2", "syn"]

    # Plot output
    fig = Figure(; resolution=(800, 200*length(stages)))
    axs = [Axis(fig[i, 1]) for i in eachindex(stages)]
    for (ax, stage) in zip(axs, stages)
        lines!(ax, resps_orig[stage]; color=:gray)
        lines!(ax, resps_new[stage]; color=:red, linestyle=:dash, linewidth=3.0)
        ax.ylabel = stage
    end

    # Adjust labels
    xlims!.(axs[1:2], 1000, 1100)
    xlims!(axs[3], 5000, 6000)
    axs[end].xlabel = "Time (s)"
    hidexdecorations!.(axs[1:(end-1)], ticks=false, grid=false)
    fig
    save(projectdir("outputs", "test", "responses_beyond_powerlaw.png"), fig)
end
plot_model_2!()
