using Helios
using AuditorySignalUtils
using AuditoryNerveFiber
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
    meout_new = zeros(length(x))
    controlout_new = zeros(length(x))
    c1out_new = zeros(length(x))
    c1vihcout_new = zeros(length(x))
    c2out_new = zeros(length(x))
    c2vihcout_new = zeros(length(x))
    ihcout_new = zeros(length(x))

    # Run model
    model!(
        x, 
        1000.0, 
        1, 
        1/100e3, 
        length(x), 
        1.0, 
        1.0, 
        2, 
        meout_new, 
        controlout_new, 
        c1out_new, 
        c1vihcout_new, 
        c2out_new, 
        c2vihcout_new, 
        ihcout_new
    )

    # Get all responses from ANF.jl debug
    ihcout = zeros(length(x))
    c1out = zeros(length(x))
    c1vihcout = zeros(length(x))
    c2out = zeros(length(x))
    c2vihcout = zeros(length(x))
    controlout = zeros(length(x))
    ccall(
        (:IHCDEBUG, "C:\\Users\\dguest2\\cl_code\\Helios\\external\\julia\\libzbc2014debug.so"),
        Cvoid,                   # return type
        (                        # arg types
            Ptr{Cdouble},        # px
            Cdouble,             # cf
            Cint,                # nrep
            Cdouble,             # tdres
            Cint,                # totalstim
            Cdouble,             # cohc
            Cdouble,             # cihc
            Cint,                # species
            Ptr{Cdouble},        # ihcout
            Ptr{Cdouble},        # c1out
            Ptr{Cdouble},        # c1vihcout
            Ptr{Cdouble},        # c1out
            Ptr{Cdouble},        # c2vihcout
            Ptr{Cdouble},        # controlout
        ),
        x, 1000.0, 1, 1/100e3, length(x), 1.0, 1.0, 2, ihcout, c1out, c1vihcout, c2out, c2vihcout, controlout, # pass arguments
    )

    # Plot output
    fig = Figure()
    axs = [Axis(fig[i, 1]) for i in 1:7]
    t = 0.0:(1/100e3):(0.25-1/100e3)
    for (ax, resp) in zip(axs, [x, controlout, c1out, c1vihcout, c2out, c2vihcout, ihcout])
        lines!(ax, t, resp; linestyle=:dash)
    end

    for (ax, resp) in zip(axs, [x, controlout_new, c1out_new, c1vihcout_new, c2out_new, c2vihcout_new, ihcout_new])
        lines!(ax, t, resp)
    end

    xlims!.(axs, 0.03, 0.04)
    axs[end].xlabel = "Time (s)"
    for (ax, label) in zip(axs, ["Input pressure", "Control", "C1", "C1 VIHC", "C2", "C2 VIHC", "IHC"])
        ax.ylabel = label
    end
    hidexdecorations!.(axs[1:(end-1)], ticks=false, grid=false)
    fig
    save(projectdir("outputs", "test", "model!.png"), fig)
end
plot_model!()
