export sim_anrate_gfc2023, sim_gfc2023, sim_gfc2023_dict, sim_orig, sim_orig_dict

"""
    sim_gfc2023(input, cf; fs=100e3, fs_synapse=10e3, fiber_type="high", power_law="approximate", fractional=false, n_rep=1)

Simulates full model output for sound-pressure input

# Arguments
- `input::Vector{Float64}`: sound-pressure waveform (Pa)
- `cf::Float64`: characteristic frequency of the fiber in Hz
- `fs::Float64`: sampling rate of the *input* in Hz
- `cohc::Float64`:
- `cihc::Float64`:
- `species::String`:
- `fiber_type::String`: fiber type, one of ("low", "medium", "high") spontaneous rate
- `fractional::Bool`: 

# Returns
- `output::Vector{Float64}`: synapse output (unknown units?), length is `length(input)`
"""
function sim_gfc2023(
    x::Vector{Float64}, 
    cf::Vector{Float64}; 
    fs::Float64=100e3,
    cohc::Float64=1.0,
    cihc::Float64=1.0,
    species::String="human",
    fiber_type::String="high", 
    fractional=false,
)
    # Calculate n_chan
    n_chan = length(cf)

    # Convert human-readable arguments into C-side floats/ints
    species_flag = Dict(
        "cat" => 1,
        "human" => 2,
        "human_glasberg" => 3
    )[species]

    spont = Dict(
        "low" => 0.1,
        "medium" => 4.0,
        "high" => 100.0
    )[fiber_type]

    # Synthesize ffGn
    if fractional
        ffGn = map(1:n_chan) do _
            ffGn_native(
                Int(ceil(length(x))),
                1/fs,
                0.9,
                1.0,
                spont,
            )
        end
    else
        ffGn = [zeros(Int(ceil(length(x)))) for _ in 1:n_chan]
    end

    # Pre-allocate memory
    controlout = [zeros(length(x)) for _ in 1:n_chan]
    c1out = [zeros(length(x)) for _ in 1:n_chan]
    c2out = [zeros(length(x)) for _ in 1:n_chan]
    ihcout = [zeros(length(x)) for _ in 1:n_chan]
    exponout = [zeros(length(x)) for _ in 1:n_chan]
    sout1 = [zeros(length(x)) for _ in 1:n_chan]
    sout2 = [zeros(length(x)) for _ in 1:n_chan]
    synout = [zeros(length(x)) for _ in 1:n_chan]
    anrateout = [zeros(length(x)) for _ in 1:n_chan]

    # Run model
    model!(
        x, 
        ffGn,
        cf,
        n_chan,
        1/fs, 
        length(x), 
        cohc, 
        cihc, 
        species_flag, 
        spont,
        controlout, 
        c1out, 
        c2out, 
        ihcout,
        exponout,
        sout1,
        sout2,
        synout,
        anrateout,
    )

    # Return
    return controlout, c1out, c2out, ihcout, exponout, sout1, sout2, synout, anrateout
end

function sim_gfc2023(x::Vector{Float64}, cf::Float64; kwargs...)
    [x[1] for (idx, x) in enumerate(sim_gfc2023(x, [cf]; kwargs...))]
end

function sim_gfc2023_dict(args...; kwargs...)
    control, c1, c2, ihc, expon, sout1, sout2, syn, anrate = sim_gfc2023(args..., kwargs...)
    return Dict(
        "control" => control,
        "c1" => c1,
        "c2" => c2,
        "ihc" => ihc,
        "expon" => expon,
        "sout1" => sout1,
        "sout2" => sout2,
        "syn" => syn,
        "anrate" => anrate,
    )
end

function sim_anrate_gfc2023(x::Vector{Float64}, cf::Float64; kwargs...)
    _, _, _, _, _, _, _, synout, _, _, _ = sim_gfc2023(x, cf; kwargs...)
    synout ./ (1.0 .+ 0.75e-3 .* synout)
end

function sim_anrate_gfc2023(x::Vector{Float64}, cf::Vector{Float64}; kwargs...)
    _, _, _, _, _, _, _, synout, _, _, _ = sim_gfc2023(x, cf; kwargs...)
    map(x -> x ./ (1.0 .+ 0.75e-3 .* x), synout)
end

function sim_orig(
    x::Vector{Float64}, 
    cf::Float64; 
    fs::Float64=100e3,
    cohc::Float64=1.0,
    cihc::Float64=1.0,
    species::String="human",
    fiber_type::String="high", 
    power_law::String="actual", 
    fractional::Bool=false,
)
    # Convert human-readable arguments into C-side floats/ints
    species_flag = Dict(
        "cat" => 1,
        "human" => 2,
        "human_glasberg" => 3
    )[species]

    spont = Dict(
        "low" => 0.1,
        "medium" => 4.0,
        "high" => 100.0
    )[fiber_type]

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
            spont,
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
        ihcout, ffGn, 1/fs, cf, length(ihcout), 1, spont, noiseType, implnt, 10e3, synout, exponout, powerlawin, sout1, sout2, @cfunction(decimate, Ptr{Cdouble}, (Ptr{Cdouble}, Cint, Cint)),
    )

    # Return
    return controlout, c1out, c1vihcout, c2out, c2vihcout, ihcout, synout, exponout, powerlawin, sout1, sout2
end

function sim_orig_dict(args...; kwargs...)
    control, c1, c1vihc, c2, c2vihc, ihc, syn, expon, powerlaw, sout1, sout2 = sim_orig(args..., kwargs...)
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
    )
end
