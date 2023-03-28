using CairoMakie
using Helios_plots
using Helios
using CarneyLabUtils2
using AuditorySignalUtils
using DSP

### ========================================================================================
# Configure range of models and stimuli we want to look at
### ========================================================================================
# Set parameters for stimulus and model
n_cf = 31
cf = LogRange(1000.0, 6000.0, n_cf)
fs = 50e3 
f0 = 100.0

# Models
weights_wdr = [0.0, 1.0, 2.0, 4.0, 8.0]
weights_ic = [0.0, 1.0, 2.0, 4.0, 8.0]
lens_integ = [1, 2, 4, 8]
models = map(Iterators.product(weights_wdr, weights_ic, lens_integ)) do (weight_wdr, weight_ic, len_integ)
    m = AuditorySubcortexGFCv1(;
        cf=cf,
        moc_weight_wdr=weight_wdr,
        moc_weight_ic=weight_ic,
        moc_len_integ=len_integ,
        moc_cutoff=0.2,
        ic_a=2.0,
        fractional=false,
        fs=fs,
        pathhead="\\home\\daniel\\cl_sim\\efferent\\explore_3-27",
    )
    return m
end

# Stimuli
stimuli = map([3000.0, 4000.0]) do f2
    s = TwoFormantTone(; f2=f2, fs=fs)
    return s
end

# Collect itr 
itr = Iterators.product(models, stimuli)

### ========================================================================================
# Look at two-formant stimulus
### ========================================================================================
for ele in itr
    # Synthesize two-formant complex tone
    m, s = ele
    if isfile(filename([m, s]))
        temp = load(filename([m, s]))
        x, resp = temp["x"], temp["resp"]
    else
        x = synthesize!(s)
        resp = m(x)
        save(filename([m, s]), Dict("x" => x, "resp" => resp))
    end

    # Plot neurograms
    neurograms = map(["ihc", "hsr", "lsr", "ic", "gain"]) do output
        fig=Figure(; resolution=(1800, 800))
        fig, ax = plot_neurogram(
            x, 
            cf, 
            resp[output], 
            fig; 
            clims=getylim(output), 
            ylims=getylim_avg(output), 
            cmap=getcmap(output)
        )
        return fig
    end
    fig = displayimg(hcat(getimg.(neurograms)...))

    # Add title to neurograms
    t1 = title(m, [:moc_weight_wdr, :moc_weight_ic, :moc_len_integ])
    t2 = title(s, [:f1, :b1, :a1, :f2, :b2, :a2])
    Label(fig[0, :], t1 * "\n" * t2; tellwidth=false)

    # Hand-construct output filename for easy filtering
    @unpack moc_weight_wdr, moc_weight_ic, moc_len_integ = m
    @unpack f2 = s
    fn1 = identifier(m, [:moc_weight_wdr, :moc_weight_ic, :moc_len_integ])
    fn2 = identifier(s, [:f2])
    fn = "twoformant_neurograms_$(fn1)_$(fn2)"
    save("\\home\\daniel\\cl_fig\\efferent\\3-27\\$fn.png", fig)
end