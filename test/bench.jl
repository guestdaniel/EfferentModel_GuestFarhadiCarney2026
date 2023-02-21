using Helios
using AuditorySignalUtils
using CairoMakie
using DrWatson

# ==========================================================================================
# Middle-ear filter
# ==========================================================================================
function test_middle_ear_filter()
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
    save(projectdir("outputs", "test", "middle_ear.png"), fig)
end
test_middle_ear_filter()