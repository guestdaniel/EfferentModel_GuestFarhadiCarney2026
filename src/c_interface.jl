export model!

# ##############################################################################
# Handle to the model library
# ##############################################################################
# The library is opened on first use and kept open, with the addresses of its entry points
# cached. Before every call, the library file's modification time is checked; if the file
# has changed (e.g., after rebuilding with deps/build_clib.sh), the library is closed and
# reopened so that the new build is used. The cache is not thread-safe.
mutable struct ModelLibrary
    handle::Ptr{Cvoid}
    mtime::Float64
    monaural::Ptr{Cvoid}
    binaural::Ptr{Cvoid}
end

const MODEL_LIBRARY = ModelLibrary(C_NULL, 0.0, C_NULL, C_NULL)

# Called from the module's __init__: pointers from a previous process are never valid
function reset_model_library!()
    MODEL_LIBRARY.handle = C_NULL
    MODEL_LIBRARY.mtime = 0.0
    MODEL_LIBRARY.monaural = C_NULL
    MODEL_LIBRARY.binaural = C_NULL
    return nothing
end

"""
    model_function(name::Symbol)

Return the address of entry point `name` (`:model_monaural` or `:model_binaural`) in the
model library, opening the library if needed and reopening it if the file has changed.
"""
function model_function(name::Symbol)
    lib = MODEL_LIBRARY
    t = mtime(libgfc2026)
    if lib.handle == C_NULL || t != lib.mtime
        lib.handle == C_NULL || Libdl.dlclose(lib.handle)
        lib.handle = C_NULL
        handle = Libdl.dlopen(libgfc2026)
        lib.monaural = Libdl.dlsym(handle, :model_monaural)
        lib.binaural = Libdl.dlsym(handle, :model_binaural)
        lib.handle = handle
        lib.mtime = t
    end
    name === :model_monaural && return lib.monaural
    name === :model_binaural && return lib.binaural
    throw(ArgumentError("unknown model entry point: $name"))
end

"""
    model!(px, ffGn_hsr, ffGn_lsr, cf, n_chan, tdres, totalstim, ..., gain, gainpostmix)

Low-level interface to the C model: runs `model_monaural` (when `px` is a
`Vector{Float64}`) or `model_binaural` (when `px` is a `Vector{Vector{Float64}}`), writing
results in place into the output arrays `controlout` through `gainpostmix`. Most users
should call `sim_gfc2023` or `sim_gfc2023!` instead.

Every output array (and every channel within it) must be a separate allocation that does
not share memory with any other argument. Do not pass the same array for several outputs
(for example, one scratch buffer for outputs that are not needed): the C code assumes the
arrays do not overlap, and overlapping arrays give silently incorrect results.
"""
function model!(
    px::Vector{Float64},
    ffGn_hsr::Vector{Vector{Float64}},
    ffGn_lsr::Vector{Vector{Float64}},
    cf::Vector{Float64},
    n_chan::Int64,
    tdres::Float64,
    totalstim::Int64, 
    cohc::Vector{Float64},
    cihc::Vector{Float64},
    species::Int64,
    spont::Float64,
    powerlaw_mode::Int64,
    cn_tau_e::Float64,
    cn_tau_i::Float64,
    cn_delay::Float64,
    cn_amp::Float64,
    cn_inh::Float64,
    ic_tau_e::Float64,
    ic_tau_i::Float64,
    ic_delay::Float64,
    ic_amp::Float64,
    ic_inh::Float64,
    moc_cutoff::Float64,
    moc_beta::Vector{Float64},
    moc_offset::Vector{Float64},
    moc_minval::Float64,
    moc_maxval::Float64,
    moc_weight::Vector{Float64},
    moc_width::Float64,
    dur_settle::Float64,
    moc_delay::Float64,
    moc_fix_gain::Int64,
    controlout::Vector{Vector{Float64}},
    c1out::Vector{Vector{Float64}},
    c2out::Vector{Vector{Float64}},
    ihcout::Vector{Vector{Float64}},
    expout_hsr::Vector{Vector{Float64}},
    sout1_hsr::Vector{Vector{Float64}},
    sout2_hsr::Vector{Vector{Float64}},
    synout_hsr::Vector{Vector{Float64}},
    expout_lsr::Vector{Vector{Float64}},
    sout1_lsr::Vector{Vector{Float64}},
    sout2_lsr::Vector{Vector{Float64}},
    synout_lsr::Vector{Vector{Float64}},
    hsrout::Vector{Vector{Float64}},
    lsrout::Vector{Vector{Float64}},
    cnout::Vector{Vector{Float64}},
    icout::Vector{Vector{Float64}},
    mocwdr::Vector{Vector{Float64}},
    mocic::Vector{Vector{Float64}},
    gain::Vector{Vector{Float64}},
    gainpostmix::Vector{Vector{Float64}},
)
    modelfunc = model_function(:model_monaural)

    # Place call
    ccall(
        modelfunc,  # pointer to function from lib
        Cvoid, # return type
        (
            Ptr{Cdouble}, # px
            Ptr{Ptr{Cdouble}}, # ffGn_hsr
            Ptr{Ptr{Cdouble}}, # ffGn_lsr
            Ptr{Cdouble}, # cf
            Cint,         # nchan
            Cdouble,      # tdres
            Cint,         # totalstim
            Ptr{Cdouble}, # cohc
            Ptr{Cdouble}, # cihc
            Cint,         # species
            Cdouble,      # spont
            Cint,         # powerlaw_mode
            Cdouble,      # cn_tau_e
            Cdouble,      # cn_tau_i,
            Cdouble,      # cn_delay
            Cdouble,      # cn_amp
            Cdouble,      # cn_inh
            Cdouble,      # ic_tau_e
            Cdouble,      # ic_tau_i,
            Cdouble,      # ic_delay
            Cdouble,      # ic_amp
            Cdouble,      # ic_inh
            Cdouble,      # moc_cutoff
            Ptr{Cdouble}, # moc_beta
            Ptr{Cdouble}, # moc_offset
            Cdouble,      # moc_minrate
            Cdouble,      # moc_maxrate
            Ptr{Cdouble}, # moc_weight
            Cdouble,      # moc_width
            Cdouble,      # dur_settle
            Cdouble,      # moc_delay
            Cint,         # moc_fix_gain
            Ptr{Ptr{Cdouble}}, # control
            Ptr{Ptr{Cdouble}}, # c1 
            Ptr{Ptr{Cdouble}}, # c2 
            Ptr{Ptr{Cdouble}}, # ihcout
            Ptr{Ptr{Cdouble}}, # expout_hsr
            Ptr{Ptr{Cdouble}}, # sout1_hsr
            Ptr{Ptr{Cdouble}}, # sout2_hsr
            Ptr{Ptr{Cdouble}}, # synout_hsr
            Ptr{Ptr{Cdouble}}, # expout_lsr
            Ptr{Ptr{Cdouble}}, # sout1_lsr
            Ptr{Ptr{Cdouble}}, # sout2_lsr
            Ptr{Ptr{Cdouble}}, # synout_lsr
            Ptr{Ptr{Cdouble}}, # hsrout
            Ptr{Ptr{Cdouble}}, # lsrout
            Ptr{Ptr{Cdouble}}, # cnout
            Ptr{Ptr{Cdouble}}, # icout
            Ptr{Ptr{Cdouble}}, # mocwdr
            Ptr{Ptr{Cdouble}}, # mocic
            Ptr{Ptr{Cdouble}}, # gain
            Ptr{Ptr{Cdouble}}, # gainpostmix
        ),
        px,
        ffGn_hsr,
        ffGn_lsr,
        cf, 
        n_chan,
        tdres, 
        totalstim, 
        cohc, 
        cihc, 
        species, 
        spont,
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
        moc_beta,
        moc_offset,
        moc_minval,
        moc_maxval,
        moc_weight,
        moc_width,
        dur_settle,
        moc_delay,
        moc_fix_gain,
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
        gainpostmix,
    )
    return nothing
end

function model!(
    px::Vector{Vector{Float64}},
    ffGn_hsr::Vector{Vector{Vector{Float64}}},
    ffGn_lsr::Vector{Vector{Vector{Float64}}},
    cf::Vector{Float64},
    n_chan::Int64,
    tdres::Float64,
    totalstim::Int64,
    cohc::Vector{Vector{Float64}},
    cihc::Vector{Vector{Float64}},
    species::Int64,
    spont::Float64,
    powerlaw_mode::Int64,
    cn_tau_e::Float64,
    cn_tau_i::Float64,
    cn_delay::Float64,
    cn_amp::Float64,
    cn_inh::Float64,
    ic_tau_e::Float64,
    ic_tau_i::Float64,
    ic_delay::Float64,
    ic_amp::Float64,
    ic_inh::Float64,
    moc_cutoff::Float64,
    moc_beta::Vector{Float64},
    moc_offset::Vector{Float64},
    moc_minval::Float64,
    moc_maxval::Float64,
    moc_weight::Vector{Float64},
    moc_width::Float64,
    dur_settle::Float64,
    moc_delay::Float64,
    moc_fix_gain::Int64,
    controlout::Vector{Vector{Vector{Float64}}},
    c1out::Vector{Vector{Vector{Float64}}},
    c2out::Vector{Vector{Vector{Float64}}},
    ihcout::Vector{Vector{Vector{Float64}}},
    expout_hsr::Vector{Vector{Vector{Float64}}},
    sout1_hsr::Vector{Vector{Vector{Float64}}},
    sout2_hsr::Vector{Vector{Vector{Float64}}},
    synout_hsr::Vector{Vector{Vector{Float64}}},
    expout_lsr::Vector{Vector{Vector{Float64}}},
    sout1_lsr::Vector{Vector{Vector{Float64}}},
    sout2_lsr::Vector{Vector{Vector{Float64}}},
    synout_lsr::Vector{Vector{Vector{Float64}}},
    hsrout::Vector{Vector{Vector{Float64}}},
    lsrout::Vector{Vector{Vector{Float64}}},
    cnout::Vector{Vector{Vector{Float64}}},
    icout::Vector{Vector{Vector{Float64}}},
    mocwdr::Vector{Vector{Vector{Float64}}},
    mocic::Vector{Vector{Vector{Float64}}},
    gain::Vector{Vector{Vector{Float64}}},
    gainpostmix::Vector{Vector{Vector{Float64}}},
)
    modelfunc = model_function(:model_binaural)

    ccall(
        modelfunc,
        Cvoid,
        (
            Ptr{Ptr{Cdouble}},          # px
            Ptr{Ptr{Ptr{Cdouble}}},     # ffGn_hsr
            Ptr{Ptr{Ptr{Cdouble}}},     # ffGn_lsr
            Ptr{Cdouble},               # cf
            Cint,                       # n_chan
            Cdouble,                    # tdres
            Cint,                       # totalstim
            Ptr{Ptr{Cdouble}},          # cohc
            Ptr{Ptr{Cdouble}},          # cihc
            Cint,                       # species
            Cdouble,                    # spont
            Cint,                       # powerlaw_mode
            Cdouble,                    # cn_tau_e
            Cdouble,                    # cn_tau_i
            Cdouble,                    # cn_delay
            Cdouble,                    # cn_amp
            Cdouble,                    # cn_inh
            Cdouble,                    # ic_tau_e
            Cdouble,                    # ic_tau_i
            Cdouble,                    # ic_delay
            Cdouble,                    # ic_amp
            Cdouble,                    # ic_inh
            Cdouble,                    # moc_cutoff
            Ptr{Cdouble},               # moc_beta
            Ptr{Cdouble},               # moc_offset
            Cdouble,                    # moc_minval
            Cdouble,                    # moc_maxval
            Ptr{Cdouble},               # moc_weight
            Cdouble,                    # moc_width
            Cdouble,                    # dur_settle
            Cdouble,                    # moc_delay
            Cint,                       # moc_fix_gain
            Ptr{Ptr{Ptr{Cdouble}}},     # controlout
            Ptr{Ptr{Ptr{Cdouble}}},     # c1out
            Ptr{Ptr{Ptr{Cdouble}}},     # c2out
            Ptr{Ptr{Ptr{Cdouble}}},     # ihcout
            Ptr{Ptr{Ptr{Cdouble}}},     # expout_hsr
            Ptr{Ptr{Ptr{Cdouble}}},     # sout1_hsr
            Ptr{Ptr{Ptr{Cdouble}}},     # sout2_hsr
            Ptr{Ptr{Ptr{Cdouble}}},     # synout_hsr
            Ptr{Ptr{Ptr{Cdouble}}},     # expout_lsr
            Ptr{Ptr{Ptr{Cdouble}}},     # sout1_lsr
            Ptr{Ptr{Ptr{Cdouble}}},     # sout2_lsr
            Ptr{Ptr{Ptr{Cdouble}}},     # synout_lsr
            Ptr{Ptr{Ptr{Cdouble}}},     # hsrout
            Ptr{Ptr{Ptr{Cdouble}}},     # lsrout
            Ptr{Ptr{Ptr{Cdouble}}},     # cnout
            Ptr{Ptr{Ptr{Cdouble}}},     # icout
            Ptr{Ptr{Ptr{Cdouble}}},     # mocwdr
            Ptr{Ptr{Ptr{Cdouble}}},     # mocic
            Ptr{Ptr{Ptr{Cdouble}}},     # gain
            Ptr{Ptr{Ptr{Cdouble}}},     # gainpostmix
        ),
        px,
        ffGn_hsr,
        ffGn_lsr,
        cf,
        n_chan,
        tdres,
        totalstim,
        cohc,
        cihc,
        species,
        spont,
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
        moc_beta,
        moc_offset,
        moc_minval,
        moc_maxval,
        moc_weight,
        moc_width,
        dur_settle,
        moc_delay,
        moc_fix_gain,
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
        gainpostmix,
    )
    return nothing
end

