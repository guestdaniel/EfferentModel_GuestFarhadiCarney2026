export sim_gfc2023ia, sim_gfc2023ia_dict

"""
    sim_gfc2023ia(input, cf; fs=100e3, fs_synapse=10e3, power_law="approximate", fractional=false, n_rep=1)

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
function sim_gfc2023ia(
    x::Vector{Float64}, 
    cf::Vector{Float64}; 
    fs::Float64=100e3,
    cohc::Float64=1.0,
    cihc::Float64=1.0,
    species::String="human",
    fractional=false,
    powerlaw_mode=1,
    cn_tau_e=0.5e-3,
    cn_tau_i=2.0e-3,
    cn_delay=1.0e-3,
    cn_amp=1.5,
    cn_inh=0.6,
    ic_tau_e=1.0e-3,
    ic_tau_i=2.0e-3,
    ic_delay=1.0e-3,
    ic_amp=4.0,
    ic_inh=0.9,
    moc_cutoff=0.64,
    moc_beta_wdr=0.01,
    moc_offset_wdr=0.0,
    moc_minrate_wdr=0.001,
    moc_maxrate_wdr=1.0,
    moc_beta_ic=0.01,
    moc_offset_ic=0.0,
    moc_minrate_ic=0.001,
    moc_maxrate_ic=1.0,
    moc_weight_wdr=0.0,
    moc_weight_ic=0.0,
    moc_width_wdr=0.0,
    dur_pad_left=0.0,
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
        ffGn_hsr = map(1:10) do _
            map(1:n_chan) do _
                ffGn_native(
                    len_total,
                    1/fs,
                    0.9,
                    1.0,
                    100.0,
                )
            end
        end
        ffGn_lsr = map(1:5) do _
            map(1:n_chan) do _
                ffGn_native(
                    len_total,
                    1/fs,
                    0.9,
                    1.0,
                    0.1,
                )
            end 
        end
    else
        ffGn_hsr = [[zeros(len_total) for _ in 1:n_chan] for _ in 1:10]
        ffGn_lsr = [[zeros(len_total) for _ in 1:n_chan] for _ in 1:5]
    end

    # Convert bool powerlaw_include_fast to integer
#    powerlaw_include_fast = Int64(powerlaw_include_fast)

    # Pre-allocate memory
    controlout = [zeros(len_total) for _ in 1:n_chan]
    c1out = [zeros(len_total) for _ in 1:n_chan]
    c2out = [zeros(len_total) for _ in 1:n_chan]
    ihcout = [zeros(len_total) for _ in 1:n_chan]
    expout_hsr = [zeros(len_total) for _ in 1:n_chan]
    sout1_hsr = [[zeros(len_total) for _ in 1:n_chan] for _ in 1:10]
    sout2_hsr = [[zeros(len_total) for _ in 1:n_chan] for _ in 1:10]
    synout_hsr = [zeros(len_total) for _ in 1:n_chan]
    expout_lsr = [zeros(len_total) for _ in 1:n_chan]
    sout1_lsr = [[zeros(len_total) for _ in 1:n_chan] for _ in 1:5]
    sout2_lsr = [[zeros(len_total) for _ in 1:n_chan] for _ in 1:5]
    synout_lsr = [zeros(len_total) for _ in 1:n_chan]
    hsrout = [zeros(len_total) for _ in 1:n_chan]
    lsrout = [zeros(len_total) for _ in 1:n_chan]
    cnout = [zeros(len_total) for _ in 1:n_chan]
    icout = [zeros(len_total) for _ in 1:n_chan]
    mocwdr = [zeros(len_total) for _ in 1:n_chan]
    mocic = [zeros(len_total) for _ in 1:n_chan]
    gain = [zeros(len_total) for _ in 1:n_chan]

    # Run model
    model_ia!(
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

function sim_gfc2023ia(x::Vector{Float64}, cf::Float64; kwargs...)
    [x[1] for (idx, x) in enumerate(sim_gfc2023ia(x, [cf]; kwargs...))]
end

function sim_gfc2023ia_dict(args...; kwargs...)
    control, c1, c2, ihc, expon, sout1, sout2, syn, hsr, lsr, cn, ic, mocwdr, mocic, gain = sim_gfc2023ia(args...; kwargs...)
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
