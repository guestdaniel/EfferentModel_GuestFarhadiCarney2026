export sim_gfc2023, sim_gfc2023_dict, sim_orig, sim_orig_dict

"""
    sim_gfc2023(input, cf; fs=100e3, fs_synapse=10e3, power_law="approximate", fractional=false, n_rep=1)

Simulates full model output for sound-pressure input

# Arguments
- `input::Vector{Float64}`: sound-pressure waveform (Pa)
- `cf::Float64`: characteristic frequency of the fiber in Hz
- `fs::Float64`: sampling rate of the *input* in Hz
- `cohc::Float64`:
- `cihc::Float64`:
- `species::String`:
- `fractional::Bool`: 

# Returns
- `output::Vector{Float64}`: synapse output (unknown units?), length is `length(input)`
"""
function sim_gfc2023(
    x::Vector{Float64}, 
    cf::Vector{Float64}; 
    fs::Float64=100e3,
    cohc::Vector{Float64}=ones(size(cf)),
    cihc::Vector{Float64}=ones(size(cf)),
    species::String="human",
    fractional=false,
    powerlaw_mode=2,
    cn_tau_e=0.5e-3,
    cn_tau_i=2.0e-3,
    cn_delay=1.0e-3,
    cn_amp=1.5,
    cn_inh=0.6,
    ic_tau_e=1.0/(10.0*64.0),  # BMF = 64 Hz
    ic_tau_i=1.0/(10.0*64.0)*1.5,
    ic_delay=1.0/(10.0*64.0)*2.0,
    ic_amp=1.0,
    ic_inh=0.9,
    moc_cutoff=0.64,
    moc_beta_wdr=0.015,
    moc_offset_wdr=10.0*4.0,
    moc_minrate_wdr=0.1,
    moc_maxrate_wdr=1.0,
    moc_beta_ic=0.015,
    moc_offset_ic=10.0*4.0,
    moc_minrate_ic=0.1,
    moc_maxrate_ic=1.0,
    moc_weight_wdr=4.0,
    moc_weight_ic=4.0,
    moc_width_wdr=0.5,
    dur_pad_left=0.2,
    moc_delay=0.025,
    clip_left=dur_pad_left == 0.0 ? false : true,
    dur_pad_right=0.0,
    clip_right=dur_pad_right == 0.0 ? false : true,
)
    # Calculate pad sizes in samples
    len_pad_left = Int(floor(dur_pad_left*fs))
    len_pad_right = Int(floor(dur_pad_right*fs))
    len_stim = length(x)
    len_total = len_pad_left + len_stim + len_pad_right

    # Pad x
    stim = vcat(zeros(len_pad_left), x, zeros(len_pad_right))

    # Calculate n_chan
    n_chan = length(cf)

    # Convert human-readable arguments into C-side floats/ints
    species_flag = Dict(
        "cat" => 1,
        "human" => 2,
        "human_glasberg" => 3
    )[species]

    # Synthesize ffGn
    if fractional
        ffGn_hsr = map(1:n_chan) do _
            ffGn_native(
                len_total,
                1/fs,
                0.9,
                1.0,
                100.0,
            )
        end
        ffGn_lsr = map(1:n_chan) do _
            ffGn_native(
                len_total,
                1/fs,
                0.9,
                1.0,
                0.1,
            )
        end
    else
        ffGn_hsr = [zeros(len_total) for _ in 1:n_chan]
        ffGn_lsr = [zeros(len_total) for _ in 1:n_chan]
    end

    # Convert bool powerlaw_include_fast to integer
#    powerlaw_include_fast = Int64(powerlaw_include_fast)

    # Pre-allocate memory
    controlout = [zeros(len_total) for _ in 1:n_chan]
    c1out = [zeros(len_total) for _ in 1:n_chan]
    c2out = [zeros(len_total) for _ in 1:n_chan]
    ihcout = [zeros(len_total) for _ in 1:n_chan]
    expout_hsr = [zeros(len_total) for _ in 1:n_chan]
    sout1_hsr = [zeros(len_total) for _ in 1:n_chan]
    sout2_hsr = [zeros(len_total) for _ in 1:n_chan]
    synout_hsr = [zeros(len_total) for _ in 1:n_chan]
    expout_lsr = [zeros(len_total) for _ in 1:n_chan]
    sout1_lsr = [zeros(len_total) for _ in 1:n_chan]
    sout2_lsr = [zeros(len_total) for _ in 1:n_chan]
    synout_lsr = [zeros(len_total) for _ in 1:n_chan]
    hsrout = [zeros(len_total) for _ in 1:n_chan]
    lsrout = [zeros(len_total) for _ in 1:n_chan]
    cnout = [zeros(len_total) for _ in 1:n_chan]
    icout = [zeros(len_total) for _ in 1:n_chan]
    mocwdr = [zeros(len_total) for _ in 1:n_chan]
    mocic = [zeros(len_total) for _ in 1:n_chan]
    gain = [zeros(len_total) for _ in 1:n_chan]

    # Run model
    model!(
        stim, 
        ffGn_hsr,
        ffGn_lsr,
        cf,
        n_chan,
        1/fs, 
        len_total, 
        cohc, 
        cihc, 
        species_flag, 
        100.0,
        powerlaw_mode,
        cn_tau_e,
        cn_tau_i,
        cn_delay,
        cn_amp,
        cn_inh,
        ic_tau_e,
        ic_tau_i,
        ic_delay,
        ic_amp,
        ic_inh,
        moc_cutoff,
        moc_beta_wdr,
        moc_offset_wdr,
        moc_minrate_wdr,
        moc_maxrate_wdr,
        moc_beta_ic,
        moc_offset_ic,
        moc_minrate_ic,
        moc_maxrate_ic,
        moc_weight_wdr,
        moc_weight_ic,
        moc_width_wdr,
        dur_pad_left,
        moc_delay,
        controlout, 
        c1out, 
        c2out, 
        ihcout,
        expout_hsr,
        sout1_hsr,
        sout2_hsr,
        synout_hsr,
        expout_lsr,
        sout1_lsr,
        sout2_lsr,
        synout_lsr,
        hsrout,
        lsrout,
        cnout,
        icout,
        mocwdr,
        mocic,
        gain,
    )

    # Return
    outputs = [controlout, c1out, c2out, ihcout, expout_hsr, sout1_hsr, sout2_hsr, synout_hsr,
               hsrout, lsrout, cnout, icout, mocwdr, mocic, gain]
    if clip_left | clip_right
        outputs = map(outputs) do output
            output = map(output) do channel
                idx_left = clip_left ? (len_pad_left+1) : 1
                idx_right = clip_right ? length(channel) - len_pad_right : length(channel)
                channel = channel[idx_left:idx_right]
            end
        end
    end
    return outputs
end

function sim_gfc2023(x::Vector{Float64}, cf::Float64; kwargs...)
    [x[1] for (idx, x) in enumerate(sim_gfc2023(x, [cf]; kwargs...))]
end

function sim_gfc2023_dict(args...; kwargs...)
    control, c1, c2, ihc, expon, sout1, sout2, syn, hsr, lsr, cn, ic, mocwdr, mocic, gain = sim_gfc2023(args...; kwargs...)
    return Dict(
        "control" => control,
        "c1" => c1,
        "c2" => c2,
        "ihc" => ihc,
        "expon" => expon,
        "sout1" => sout1,
        "sout2" => sout2,
        "syn" => syn,
        "hsr" => hsr,
        "lsr" => lsr,
        "cn" => cn,
        "ic" => ic,
        "mocwdr" => mocwdr,
        "mocic" => mocic,
        "gain" => gain,
    )
end

function sim_orig(
    x::Vector{Float64}, 
    cf::Float64; 
    fs::Float64=100e3,
    cohc::Float64=1.0,
    cihc::Float64=1.0,
    species::String="human",
    power_law::String="actual", 
    fractional::Bool=false,
)
    # Convert human-readable arguments into C-side floats/ints
    species_flag = Dict(
        "cat" => 1,
        "human" => 2,
        "human_glasberg" => 3
    )[species]

    implnt = Dict(
        "actual" => 1.0,
        "approximate" => 0.0
    )[power_law]

    noiseType = Dict(
        true => 1.0,
        false => 0.0
    )[fractional]

    # Synthesize ffGn
    if noiseType == 1.0
        ffGn = ffGn_native(
            Int(ceil((length(x) + 2 * floor(7500 / (cf / 1e3))) * 1/fs * 10e3)),
            1/fs_synapse,
            0.9,
            noiseType,
            100.0,
        )
    else
        ffGn = zeros(Int(ceil((length(x) + 2 * floor(7500 / (cf / 1e3))) * 1/fs * 10e3)))
    end

    # Call original model functions
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
        x, cf, 1, 1/fs, length(x), cohc, cihc, species_flag, ihcout, c1out, c1vihcout, c2out, c2vihcout, controlout, # pass arguments
    )

    synout = zeros(length(ihcout))
    exponout = zeros(length(x))
    delaypoint = Int(floor(7500 / (cf / 1e3)))
    powerlawin = zeros(length(x) + delaypoint*3)
    sout1 = zeros(Int(ceil((length(ihcout)+2*delaypoint) * 1/100e3 * 10e3)))
    sout2 = zeros(Int(ceil((length(ihcout)+2*delaypoint) * 1/100e3 * 10e3)))
    len_noise = Int(ceil((length(ihcout) + 2 * floor(7500 / (cf / 1e3))) * 1/fs * 10e3))
    ffGn = zeros(len_noise)
    ccall(
        (:SYNAPSEDEBUG, "C:\\Users\\dguest2\\cl_code\\Helios\\external\\julia\\libzbc2014debug.so"),
        Cvoid,                   # return type
        (                        # arg types
            Ptr{Cdouble},        # px
            Ptr{Cdouble},        # randNums
            Cdouble,             # tdres
            Cdouble,             # cf
            Cint,                # totalstim
            Cint,                # nrep
            Cdouble,             # spont
            Cdouble,             # noisetype
            Cdouble,             # implementation
            Cdouble,             # sampFreq
            Ptr{Cdouble},        # synout
            Ptr{Cdouble},        # exponout
            Ptr{Cdouble},        # powerlawin
            Ptr{Cdouble},        # sout1
            Ptr{Cdouble},        # sout2
            Ptr{Cvoid},          # decimate function handle
        ),
        ihcout, ffGn, 1/fs, cf, length(ihcout), 1, 100.0, noiseType, implnt, 10e3, synout, exponout, powerlawin, sout1, sout2, @cfunction(decimate, Ptr{Cdouble}, (Ptr{Cdouble}, Cint, Cint)),
    )

    if noiseType == 1.0
        ffGn = ffGn_native(
            Int(ceil((length(x) + 2 * floor(7500 / (cf / 1e3))) * 1/fs * 10e3)),
            1/fs_synapse,
            0.9,
            noiseType,
            0.1,
        )
    else
        ffGn = zeros(Int(ceil((length(x) + 2 * floor(7500 / (cf / 1e3))) * 1/fs * 10e3)))
    end
    synout_lsr = zeros(length(ihcout))
    exponout_lsr = zeros(length(x))
    powerlawin_lsr = zeros(length(x) + delaypoint*3)
    sout1_lsr = zeros(Int(ceil((length(ihcout)+2*delaypoint) * 1/100e3 * 10e3)))
    sout2_lsr = zeros(Int(ceil((length(ihcout)+2*delaypoint) * 1/100e3 * 10e3)))
    ccall(
        (:SYNAPSEDEBUG, "C:\\Users\\dguest2\\cl_code\\Helios\\external\\julia\\libzbc2014debug.so"),
        Cvoid,                   # return type
        (                        # arg types
            Ptr{Cdouble},        # px
            Ptr{Cdouble},        # randNums
            Cdouble,             # tdres
            Cdouble,             # cf
            Cint,                # totalstim
            Cint,                # nrep
            Cdouble,             # spont
            Cdouble,             # noisetype
            Cdouble,             # implementation
            Cdouble,             # sampFreq
            Ptr{Cdouble},        # synout
            Ptr{Cdouble},        # exponout
            Ptr{Cdouble},        # powerlawin
            Ptr{Cdouble},        # sout1
            Ptr{Cdouble},        # sout2
            Ptr{Cvoid},          # decimate function handle
        ),
        ihcout, ffGn, 1/fs, cf, length(ihcout), 1, 0.1, noiseType, implnt, 10e3, synout_lsr, exponout_lsr, powerlawin_lsr, sout1_lsr, sout2_lsr, @cfunction(decimate, Ptr{Cdouble}, (Ptr{Cdouble}, Cint, Cint)),
    )

    hsr = synout ./ (1.0 .+ 0.75e-3 .* synout)
    lsr = synout_lsr ./ (1.0 .+ 0.75e-3 .* synout_lsr)

    # Return
    return controlout, c1out, c1vihcout, c2out, c2vihcout, ihcout, synout, exponout, powerlawin, sout1, sout2, hsr, lsr
end

function sim_orig_dict(args...; kwargs...)
    control, c1, c1vihc, c2, c2vihc, ihc, syn, expon, powerlaw, sout1, sout2, hsr, lsr = sim_orig(args..., kwargs...)
    return Dict(
        "control" => control,
        "c1" => c1,
        "c1vihc" => c1vihc,
        "c2" => c2,
        "c2vihc" => c2vihc,
        "ihc" => ihc,
        "syn" => syn,
        "expon" => expon,
        "powerlaw" => powerlaw,
        "sout1" => sout1,
        "sout2" => sout2,
        "hsr" => hsr,
        "lsr" => lsr,
    )
end
