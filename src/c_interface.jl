export model!, model_wrapper!

function model_wrapper!(
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
    fastmode::Int64,
    ic_tau_e::Float64,
    ic_tau_i::Float64,
    ic_delay::Float64,
    ic_amp::Float64,
    ic_inh::Float64,
    moc_beta_wdr::Float64,
    moc_offset_wdr::Float64,
    moc_beta_ic::Float64,
    moc_offset_ic::Float64,
    moc_weight_wdr::Float64,
    moc_weight_ic::Float64,
    moc_width_wdr::Float64,
    ihcout::Vector{Vector{Float64}},
    hsrout::Vector{Vector{Float64}},
    lsrout::Vector{Vector{Float64}},
    icout::Vector{Vector{Float64}},
    gain::Vector{Vector{Float64}},
)
    ccall(
            (:model_efferent_wrapper, "C:\\Users\\dguest2\\cl_code\\Helios\\src\\model\\libgfc2023.so"), 
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
                Cint,         # fastmode
                Cdouble,      # ic_tau_e
                Cdouble,      # ic_tau_i,
                Cdouble,      # ic_delay
                Cdouble,      # ic_amp
                Cdouble,      # ic_inh
                Cdouble,      # moc_beta
                Cdouble,      # moc_offset
                Cdouble,      # moc_beta
                Cdouble,      # moc_offset
                Cdouble,      # moc_weight_wdr
                Cdouble,      # moc_weight_ic
                Cdouble,      # moc_width_wdr
                Ptr{Ptr{Cdouble}}, # ihcout
                Ptr{Ptr{Cdouble}}, # hsrout
                Ptr{Ptr{Cdouble}}, # lsrout
                Ptr{Ptr{Cdouble}}, # icout
                Ptr{Ptr{Cdouble}}, # gain
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
            fastmode,
            ic_tau_e,
            ic_tau_i,
            ic_delay,
            ic_amp,
            ic_inh,
            moc_beta_wdr,
            moc_offset_wdr,
            moc_beta_ic,
            moc_offset_ic,
            moc_weight_wdr,
            moc_weight_ic,
            moc_width_wdr,
            ihcout,
            hsrout,
            lsrout,
            icout,
            gain,
        )
end


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
    moc_beta_wdr::Float64,
    moc_offset_wdr::Float64,
    moc_minrate_wdr::Float64,
    moc_maxrate_wdr::Float64,
    moc_beta_ic::Float64,
    moc_offset_ic::Float64,
    moc_minrate_ic::Float64,
    moc_maxrate_ic::Float64,
    moc_weight_wdr::Float64,
    moc_weight_ic::Float64,
    moc_width_wdr::Float64,
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
)
    # Open library using Libdl (is there any overhead here?)
    lib = Libdl.dlopen(joinpath("src", "model", "libgfc2023.so"))
    modelfunc = Libdl.dlsym(lib, :model)

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
            Cdouble,      # moc_beta
            Cdouble,      # moc_offset
            Cdouble,      # moc_minrate
            Cdouble,      # moc_maxrate
            Cdouble,      # moc_beta
            Cdouble,      # moc_offset
            Cdouble,      # moc_minrate
            Cdouble,      # moc_maxrate
            Cdouble,      # moc_weight_wdr
            Cdouble,      # moc_weight_ic
            Cdouble,      # moc_width_wdr
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

    # Close library so it can be reloaded
    Libdl.dlclose(lib)
end

