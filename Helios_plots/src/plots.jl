export plot_timecourse!, plot_timecourse_stack!

function setylim(output::String)
    if     output == "hsr"
        return 1000.0
    elseif output == "lsr"
        return 100.0
    elseif output == "ic"
        return 1250.0
    end
end

function plot_timecourse!(
    resp::Dict, 
    channel::Int, 
    output::String="hsr", 
    fig=Figure(); 
    fs=100e3,
    title=true,
    axislabels=true,
    autoylim=false,
)
    # Derive duration and time axis
    dur = length(resp[output][channel])/fs
    t = 0.0:(1/fs):nextfloat(dur-1/fs)
    dur_moc = length(resp["gain"][channel])/fs
    t_moc = 0.0:(1/fs):nextfloat(dur_moc-1/fs)

    # Create plot
    ax = Axis(fig[1, 1])
    ax_r = Axis(fig[1, 1]; yaxisposition=:right)

    # Adjust spines, spacing, grid for overlay axis
    hidespines!(ax_r)
    hidexdecorations!(ax_r)
    hideydecorations!(ax_r, ticks=false, ticklabels=false, label=false)

    # Adjust ylimits
    ylims!(ax_r, 0.0, 1.02)
    if autoylim
        ylims!.(ax, 0.0, 1.2 * maximum(resp[output][channel]))
    else
        ylims!.(ax, 0.0, setylim(output))
    end

    # Adjust xlimits
    xlims!.([ax, ax_r], 0.0, maximum(length.([resp[k][channel] for k in keys(resp)])./fs))

    # Plot inputs
    lines!(ax, t, resp[output][channel]; color=:black)

    # Plot gain
    lines!(ax_r, t_moc, resp["gain"][channel]; color=:pink, linewidth=3.0)

    # Adjust labels and finer details
    if axislabels
        ax.xlabel = "Time re: simulation start (s)"
        ax.ylabel = "Output signal"
        ax_r.ylabel = "Cochlear \"gain\" (proportion)"
    end
    if title
        ax.title = output
    end
    fig
end

function plot_timecourse!(resp::Dict, output::Vector{String}=["hsr"], fig=Figure(); kwargs...)
    # Loop through outputs and create
    for (idx, o) in enumerate(output)
        plot_timecourse!(resp, o, fig[idx, 1]; kwargs...)
    end

    # Return fig
    fig
end

function plot_timecourse!(resp::Dict, channel::Int, output::Vector{String}=["hsr"], fig=Figure(); kwargs...)
    # Loop through outputs and create
    for (idx, o) in enumerate(output)
        plot_timecourse!(resp, channel, o, fig[idx, 1]; kwargs...)
    end

    # Return fig
    fig
end

function plot_timecourse_stack!(
    cf::Vector,
    resp::Dict, 
    output::Vector{String}=["hsr"], 
    channels=1:length(resp[output[1]]),
    fig=Figure();
    kwargs...
)
    # Loop through outputs and create
    for (idx_chan, chan) in enumerate(reverse(channels))
        for (idx_output, o) in enumerate(output)
            plot_timecourse!(resp, chan, o, fig[idx_chan, idx_output]; title=false, axislabels=false, kwargs...)
        end
    end

    # Loop through channels and add labels
    for (idx_chan, chan) in enumerate(reverse(channels))
        Label(fig[idx_chan, 0], "Channel #$chan\nCF=$(round(cf[chan])) Hz"; tellheight=false)
    end

    resize_to_layout!(fig)

    # Return fig
    fig
end
