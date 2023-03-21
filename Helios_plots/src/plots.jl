export plot_timecourse!

"""
    plot_timecourse!(resp::Dict, output::String, fig, ax)

# Arguments
- `resp::Dict`: Dict containing model response from `sim_gfc2023_dict`
- `output::String`: String indicating which stage to visualize
"""
function plot_timecourse!(resp::Dict, output::String="hsr", fig=Figure(); fs=100e3)
    # Derive duration and time axis
    dur = length(resp[output])/fs
    t = 0.0:(1/fs):nextfloat(dur-1/fs)
    dur_moc = length(resp["moc"])/fs
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
    ylims!.(ax, 0.0, 1.2 * maximum(resp[output]))

    # Adjust xlimits
    xlims!.([ax, ax_r], 0.0, maximum(length.([resp[k] for k in keys(resp)])./fs))

    # Plot inputs
    lines!(ax, t, resp[output]; color=:black)

    # Plot gain
    lines!(ax_r, t_moc, resp["moc"]; color=:pink, linewidth=3.0)

    # Adjust labels and finer details
    ax.xlabel = "Time re: simulation start (s)"
    ax.ylabel = "Output signal"
    ax_r.ylabel = "Cochlear \"gain\" (proportion)"
    ax.title = output
    fig
end

"""
    plot_timecourse!(resp::Dict, output::Vector{String}, fig, ax)

# Arguments
- `resp::Dict`: Dict containing model response from `sim_gfc2023_dict`
- `output::Vector{String}`: Vector of strings indicating which stage(s) to visualize
"""
function plot_timecourse!(resp::Dict, output::Vector{String}=["hsr"], fig=Figure(); kwargs...)
    # Loop through outputs and create
    for (idx, o) in enumerate(output)
        plot_timecourse!(resp, o, fig[idx, 1]; kwargs...)
    end

    # Return fig
    fig
end
