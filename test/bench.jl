using Helios
using AuditorySignalUtils
using CairoMakie
using DrWatson

# ==========================================================================================
# Middle-ear filter
# ==========================================================================================
function plot_middle_ear!()
    # Create inputs
    x = scale_dbspl(pure_tone(1000.0, 0.0, 0.25, 100e3), 50.0)
    meout = zeros(length(x))

    # Apply ME filter
    middle_ear!(x, 1/100e3, length(x), 1, meout)

    # Plot output
    fig = Figure()
    ax = Axis(fig[1, 1])
    t = 0.0:(1/100e3):(0.25-1/100e3)
    lines!(ax, t, x)
    lines!(ax, t, meout)
    xlims!(ax, 0.00, 0.01)
    ax.xlabel = "Time (s)"
    ax.ylabel = "Pressure (Pa)?"
    fig
    save(projectdir("outputs", "test", "middle_ear!.png"), fig)
end
plot_middle_ear!()

# ==========================================================================================
# model!
# ==========================================================================================
function plot_model!()
    # Create inputs
    x = scale_dbspl(pure_tone(1000.0, 0.0, 0.25, 100e3), 50.0)
    meout = zeros(length(x))

    # Apply ME filter
    model!(x, 1000.0, 1, 1/100e3, length(x), 1.0, 1.0, 1, meout)

    # Plot output
    fig = Figure()
    axs = [Axis(fig[i, 1]) for i in 1:2]
    t = 0.0:(1/100e3):(0.25-1/100e3)
    lines!(axs[1], t, x)
    lines!(axs[2], t, meout)
    xlims!.(axs, 0.00, 0.01)
    axs[2].xlabel = "Time (s)"
    axs[1].ylabel = "Input pressure"
    axs[2].ylabel = "Middle ear response"
    fig
    save(projectdir("outputs", "test", "model!.png"), fig)
end
plot_model!()