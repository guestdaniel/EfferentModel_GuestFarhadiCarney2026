using CairoMakie
using Helios_plots
using Helios
using CarneyLabUtils2
using AuditorySignalUtils
using DSP

### ========================================================================================
### Write code for setup and processing
### ========================================================================================

function setup(;
    moc_weight_wdr=1.0,
    moc_weight_ic=1.0,
    moc_len_integ=1,
    cf_low=1000.0,
    cf_high=6000.0,
    n_chan_per_oct=10,
    fs=50e3,
    kwargs...
)
    # Select CFs
    n_cf = Int(round((log2(cf_high) - log2(cf_low)) * n_chan_per_oct))
    cf = LogRange(cf_low, cf_high, n_cf)

    # Set up model and stimulus
    m = AuditorySubcortexGFCv1(;
        cf=cf,
        moc_weight_wdr=moc_weight_wdr,
        moc_weight_ic=moc_weight_ic,
        moc_len_integ=moc_len_integ,
        fs=fs,
        pathhead="\\home\\daniel\\cl_sim\\efferent\\3-28",
        kwargs...
    )
    s = TwoFormantTone(; fs=fs)
    return s, m
end

function run_sims(
    stages=["ihc", "hsr", "lsr", "ic", "gain"],
    resolution=(1500, 800);
    kwargs...
)
    # Setup and simulate responses
    s, m = setup(; kwargs...)

    # Load from cache if available
    if isfile(filename([m, s]))
        temp = load(filename([m, s]))
        x, r = temp["x"], temp["r"]
    else
        x = synthesize!(s)
        r = m(x)
        save(filename([m, s]), Dict("x" => x, "r" => r))
    end

    # Plot neurograms
    neurograms = map(stages) do output
        fig=Figure(; resolution=resolution)
        fig, ax = plot_neurogram(
            x, 
            m.cf, 
            r[output], 
            fig; 
            clims=getylim(output), 
            ylims=getylim_avg(output), 
            cmap=getcmap(output),
            fs=fs,
        )
        return fig
    end
    fig = displayimg(hcat(getimg.(neurograms)...))

    # Add title to neurograms
    t1 = title(m, [:moc_weight_wdr, :moc_weight_ic, :moc_len_integ])
    t2 = title(s, [:f1, :b1, :a1, :f2, :b2, :a2])
    Label(fig[0, :], t1 * "\n" * t2; tellwidth=false, fontsize=40.0)
    fig
end

### ========================================================================================
### Figure 1: Systematically vary moc_len_integ @ 10 CFs per octave
### ========================================================================================
figs = map(x -> run_sims(; moc_len_integ=x, moc_weight_ic=0.0), [0, 1, 2, 4, 8, 16])
fig = displayimg(vcat(getimg.(figs)...))
save("\\home\\daniel\\cl_fig\\efferent\\3-28\\effect_of_moc_len_integ.png", fig)

### ========================================================================================
### Figure 2: Systematically vary MOCIC strength with moc_len_integ=8
### ========================================================================================
figs = map(x -> run_sims(; moc_len_integ=8, moc_weight_ic=x), [0.0, 0.5, 1.0, 1.5, 2.0, 2.5])
fig = displayimg(vcat(getimg.(figs)...))
save("\\home\\daniel\\cl_fig\\efferent\\3-28\\effect_of_ic_strength.png", fig)

### ========================================================================================
### Figure 3: Systematically vary MOC nonlinearity parameters moc_len_integ=8, moc_weight_ic=1.0
### ========================================================================================
# Generate stackplot
figs = map(x -> run_sims(["hsr", "lsr", "mocwdr", "ic", "mocic", "gain"]; moc_len_integ=8, moc_weight_ic=1.0, moc_offset_ic=x), [0.0, 12.5, 25.0, 37.5, 50.0])
fig = displayimg(vcat(getimg.(figs)...))
save("\\home\\daniel\\cl_fig\\efferent\\3-28\\effect_of_moc_offset_ic.png", fig)

# Look a bit at MOC nonlinearity to help pick parameters
fig = Figure()
ax_wdr = Axis(fig[1, 1])
ax_ic = Axis(fig[1, 2])
x = LinRange(0.0, 200.0, 1000)
lines!(ax_wdr, x, moc_nonlinearity.(x, 0.05, 0.0, 1.0, 0.0))
lines!(ax_ic, x, moc_nonlinearity.(x, 0.2, 25.0, 1.0, 0.0))
fig

### ========================================================================================
### Figure 3: Try to hand-tune a few nice examples
### ========================================================================================
# Generate stackplot
fig = run_sims(
    ["ic", "gain"],
    (1400, 800); 
    moc_len_integ=8, 
    moc_weight_ic=1.0, 
    moc_beta_wdr=0.01,
    moc_offset_ic=35.0,
    moc_beta_ic=0.03,
)
