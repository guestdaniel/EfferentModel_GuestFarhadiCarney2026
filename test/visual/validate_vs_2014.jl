using Helios
using AuditorySignalUtils
using AuditoryNerveFiber
using DrWatson
using CairoMakie
using DSP

pt(f=1000.0, l=50.0, dur=0.2, fs=100e3) = scale_dbspl(pure_tone(f, 0.0, dur, fs), l)
update_theme!(fontsize=70)

# Create function to plot models easily
function plot_model!(name; xlims=(5000, 7000))
    # Stimulus
    x = pt(1000.0, 50.0)

    # Responses
    resps_orig = sim_orig_dict(x, 1000.0)
    resps_new = sim_gfc2023_dict(x, 1000.0)
    stages = ["control", "c1", "c2", "ihc", "expon", "sout1", "sout2", "syn"]
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
    fig = Figure(; resolution=(900*2, 270*length(stages)))
    axs = [Axis(fig[i, 1]) for i in eachindex(stages)]
    axs_diff = [Axis(fig[i, 2]) for i in eachindex(stages)]
    for (ax, ax_diff, stage) in zip(axs, axs_diff, stages)
        # Extract original and new responses
        orig = resps_orig[stage]
        new = resps_new[stage]

        # Shift by transmission delay (if needed)
        if stage in ["control", "c1", "c2"]
            orig = shiftsignal(orig, delaypoint)
        end
        
        # Resample and adjust indexing (if needed)
        if stage in ["sout1", "sout2"]
            orig = resample(orig, 10)
            delaypoint = Int(floor(7500/(1000.0/1e3)))
            orig = orig[(1:length(x)) .+ delaypoint]
        end

        lines!(ax, orig; color=:cyan, linewidth=8.0)
        lines!(ax, new; color=:red, linestyle=:dot, linewidth=3.0)
        lines!(ax_diff, new .- orig; color=:black)
        ax.ylabel = stage
    end

    # Adjust labels
    xlims!.(axs, xlims...)
    xlims!.(axs_diff, xlims...)
    axs[end].xlabel = "Samples"
    axs_diff[end].xlabel = "Samples"
    hidexdecorations!.(axs[1:(end-1)], ticks=false, grid=false)
    hidexdecorations!.(axs_diff[1:(end-1)], ticks=false, grid=false)
    save(projectdir("outputs", "test", "$name.png"), fig)
    fig
end
plot_model!("validate_sustained_responses"; xlims=(5000, 5500))
plot_model!("validate_onset_responses"; xlims=(0, 2000))

# Create function to plot models easily
function plot_model_multichannel!(name; xlims=(5000, 7000))
    # Stimulus
    x = pt(1000.0, 50.0)

    # Responses
    stages = ["control", "c1", "c2", "ihc", "expon", "sout1", "sout2", "syn"]
    resps_orig = map(_cf -> sim_orig_dict(x, _cf), [1000.0, 2000.0])
    resps_new = sim_gfc2023_dict(x, [1000.0, 2000.0])
    resps_orig = Dict([stage => map(idx -> resps_orig[idx][stage], 1:length([1000.0, 2000.0])) for stage in stages]...)

    # Plot output
    fig = Figure(; resolution=(2400*2, 600*length(stages)*2))
    axs = [Axis(fig[i, j]) for i in eachindex(stages), j in 1:2]
    for (idx_stage, stage) in enumerate(stages)
        for (idx_cf, cf) in enumerate([1000.0, 2000.0])
            delay = ccall(
                (:delay_cat, "C:\\Users\\dguest2\\cl_code\\Helios\\src\\model\\libgfc2023.so"),
                Cdouble,
                (
                    Cdouble,
                ),
                cf,
            )
            delaypoint = Int(ceil(delay/(1/100e3)))

            # Extract original and new responses
            orig = resps_orig[stage][idx_cf]
            new = resps_new[stage][idx_cf]

            # Shift by transmission delay (if needed)
            if stage in ["control", "c1", "c2"]
                orig = shiftsignal(orig, delaypoint)
            end
            
            # Resample and adjust indexing (if needed)
            if stage in ["sout1", "sout2"]
                orig = resample(orig, 10)
                delaypoint = Int(floor(7500/(cf/1e3)))
                orig = orig[(1:length(x)) .+ delaypoint]
            end

            lines!(axs[idx_stage, idx_cf], orig; color=:gray)
            lines!(axs[idx_stage, idx_cf], new; color=:red, linestyle=:dot, linewidth=2.0)
            axs[idx_stage, idx_cf].ylabel = stage
        end
    end

    # Adjust labels
    xlims!.(axs, xlims...)
    #hidexdecorations!.(axs[1:(end-1)], ticks=false, grid=false)
    save(projectdir("outputs", "test", "$name.png"), fig)
    fig
end
plot_model_multichannel!("validate_multichannel");

function plot_comparison(x1, x2; xlims=(0, length(x1)))
    fig = Figure()
    ax = Axis(fig[1, 1])
    lines!(ax, x1)
    lines!(ax, x2)
    xlims!(ax, xlims...)
    fig
end