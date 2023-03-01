export sim_gfc2023, sim_gfc2023_dict, sim_orig, sim_orig_dict

"""
    sim_gfc2023(input, cf; fs=100e3, fs_synapse=10e3, fiber_type="high", power_law="approximate", fractional=false, n_rep=1)

Simulates full model output for sound-pressure input

# Arguments
- `input::Vector{Float64}`: sound-pressure waveform (Pa)
- `cf::Float64`: characteristic frequency of the fiber in Hz
- `fs::Float64`: sampling rate of the *input* in Hz
- `fs_synapse::Float64`: sampling rate of the interior synapse simulation. The ratio between fs and fs_synapse must be an integer.
- `fiber_type::String`: fiber type, one of ("low", "medium", "high") spontaneous rate
- `power_law::String`: whether we use true or approximate power law adaptation, one of ("actual", "approximate")
- `fractional::Bool`: whether we use ffGn or not, one of (true, false)
- `n_rep::Int64`: number of repetititons to run. We assume that the input was also generated using `n_rep=n_rep`, hence we infer that the input acoustic waveform is of length `length(input)/n_rep`.

# Returns
- `output::Vector{Float64}`: synapse output (unknown units?), length is `length(input)`
"""
function sim_gfc2023(
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
            Int(ceil(length(x) + 2 * floor(7500 / (cf / 1e3)))),
            1/fs,
            0.9,
            noiseType,
            spont,
        )
    else
        ffGn = zeros(Int(ceil(length(x) + 2 * floor(7500 / (cf / 1e3)))))
    end

    # Pre-allocate memory
    meout = zeros(length(x))
    controlout = zeros(length(x))
    c1out = zeros(length(x))
    c1vihcout = zeros(length(x))
    c2out = zeros(length(x))
    c2vihcout = zeros(length(x))
    ihcout = zeros(length(x))
    synout = zeros(length(x))
    exponout = zeros(length(x))
    delaypoint = Int(floor(7500 / (cf / 1e3)))
    powerlawin = zeros(length(x) + delaypoint*3)
    sout1 = zeros(length(ihcout) + 2*delaypoint)
    sout2 = zeros(length(ihcout) + 2*delaypoint)

    # Run model
    model!(
        x, 
        ffGn,
        1000.0, 
        1/fs, 
        length(x), 
        cohc, 
        cihc, 
        species_flag, 
        spont,
        noiseType,
        implnt,
        meout, 
        controlout, 
        c1out, 
        c1vihcout, 
        c2out, 
        c2vihcout, 
        ihcout,
        synout,
        exponout,
        powerlawin,
        sout1,
        sout2,
    )

    # Return
    return meout, controlout, c1out, c1vihcout, c2out, c2vihcout, ihcout, synout, exponout, powerlawin, sout1, sout2
end

function sim_gfc2023_dict(args...; kwargs...)
    me, control, c1, c1vihc, c2, c2vihc, ihc, syn, expon, powerlaw, sout1, sout2 = sim_gfc2023(args..., kwargs...)
    return Dict(
        "me" => me,
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
        x, 1000.0, 1, 1/100e3, length(x), cohc, cihc, species_flag, ihcout, c1out, c1vihcout, c2out, c2vihcout, controlout, # pass arguments
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
        ihcout, ffGn, 1/100e3, 1000.0, length(ihcout), 1, spont, noiseType, implnt, 10e3, synout, exponout, powerlawin, sout1, sout2, @cfunction(decimate, Ptr{Cdouble}, (Ptr{Cdouble}, Cint, Cint)),
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
