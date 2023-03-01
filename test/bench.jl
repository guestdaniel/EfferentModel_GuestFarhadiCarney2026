using Helios
using AuditorySignalUtils
using AuditoryNerveFiber
using DrWatson
using CairoMakie
using DSP

pt(f=1000.0, l=50.0, dur=0.2, fs=100e3) = scale_dbspl(pure_tone(f, 0.0, dur, fs), l)

function plot_model!()
    # Stimulus
    x = pt(1000.0, 50.0)

    # Responses
    resps_orig = sim_orig_dict(x, 1000.0)
    resps_new = sim_gfc2023_dict(x, 1000.0)
    stages = ["control", "c1", "c2", "ihc", "expon"]
    delay = ccall(
        (:delay_cat, "C:\\Users\\dguest2\\cl_code\\Helios\\src\\model\\libgfc2023.so"),
        Cdouble,
        (
            Cdouble,
        ),
        1000.0
    )
    delaypoint = Int(ceil(delay/(1/100e3)))

    # Plot output
    fig = Figure(; resolution=(800, 200*length(stages)))
    axs = [Axis(fig[i, 1]) for i in eachindex(stages)]
    for (ax, stage) in zip(axs, stages)
        # Extract original and new responses
        orig = resps_orig[stage]
        new = resps_new[stage]

        # Shift by delaypoint (if needed)
        if stage in ["control", "c1", "c2"]
            orig = shiftsignal(orig, delaypoint)
        end
        lines!(ax, orig; color=:gray)
        lines!(ax, new; color=:red, linestyle=:dash, linewidth=3.0)
        ax.ylabel = stage
    end

    # Adjust labels
    xlims!.(axs, 5000, 6000)
    axs[end].xlabel = "Time (s)"
    hidexdecorations!.(axs[1:(end-1)], ticks=false, grid=false)
    fig
    save(projectdir("outputs", "test", "responses_up_to_exponout.png"), fig)
end
plot_model!()

function plot_model_2!()
    # Stimulus
    x = pt(1000.0, 50.0)

    # Responses
    resps_orig = sim_orig_dict(x, 1000.0)
    resps_new = sim_gfc2023_dict(x, 1000.0)
    stages = ["powerlaw", "sout1", "sout2", "syn"]

    # Plot output
    fig = Figure(; resolution=(800, 200*length(stages)))
    axs = [Axis(fig[i, 1]) for i in eachindex(stages)]
    for (ax, stage) in zip(axs, stages)
        lines!(ax, resps_orig[stage]; color=:gray)
        lines!(ax, resps_new[stage]; color=:red, linestyle=:dash, linewidth=3.0)
        ax.ylabel = stage
    end

    # Adjust labels
    xlims!(axs[1], 10000, 11000)
    xlims!.(axs[2:3], 700, 1000)
    xlims!(axs[4], 5000, 6000)
    axs[end].xlabel = "Time (s)"
    fig
    save(projectdir("outputs", "test", "responses_beyond_exponout.png"), fig)
end
plot_model_2!()
