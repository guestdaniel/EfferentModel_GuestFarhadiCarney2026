export AuditorySubcortexGFCv1, R4, TwoFormantTone, synthesize!, FullResponse

@with_kw struct AuditorySubcortexGFCv1 <: CarneyLabUtils2.Model
    # Params interface
    name::String="subcortex_gfc_v1"
    gitcommit::String=DrWatson.gitdescribe(projectdir())
    gitcommit_specific::Bool=false
    pathhead::String="home/daniel/cl_sim"
    exclude_identifier_keys::Vector{Symbol}=Symbol[]
    display_keys::Vector{Symbol}=[:name, :species, :fractional, :cf_low, :cf_high, :n_chan]

    # Model interface
    fs::Float64=100e3
    cf::Vector{Float64}=[1000.0]
    n_chan::Int64=length(cf)
    cf_low::Float64=minimum(cf)
    cf_high::Float64=maximum(cf)
    stage::String="ic"

    # IHC parameters
    species::String="human"

    # ANF parameters
    fractional::Bool=false

    # Subcortical parameters
    cn_τₑ::Float64=0.5e-3
    cn_τᵢ::Float64=2.0e-3
    cn_d::Float64=1.0e-3
    cn_a::Float64=1.5
    cn_s::Float64=0.6
    ic_τₑ::Float64=1.0e-3
    ic_τᵢ::Float64=2.0e-3
    ic_d::Float64=1.0e-3
    ic_a::Float64=2.0
    ic_s::Float64=0.9
    
    # Efferent parameters
    moc_cutoff::Float64=0.2
    moc_beta_wdr::Float64=0.01
    moc_offset_wdr::Float64=0.0
    moc_maxrate_wdr::Float64=1.0
    moc_minrate_wdr::Float64=0.0
    moc_beta_ic::Float64=0.01
    moc_offset_ic::Float64=0.0
    moc_maxrate_ic::Float64=1.0
    moc_minrate_ic::Float64=0.0
    moc_weight_wdr::Float64=0.0
    moc_weight_ic::Float64=0.0
    moc_len_integ::Int64=2
end

function (m::AuditorySubcortexGFCv1)(x::AbstractVector{Float64})
    sim_gfc2023_dict(
        x, 
        m.cf;
        species=m.species,
        fractional=m.fractional,
        cn_tau_e=m.cn_τₑ,
        cn_tau_i=m.cn_τᵢ,
        cn_delay=m.cn_d,
        cn_amp=m.cn_a,
        cn_inh=m.cn_s,
        ic_tau_e=m.ic_τₑ,
        ic_tau_i=m.ic_τᵢ,
        ic_delay=m.ic_d,
        ic_amp=m.ic_a,
        ic_inh=m.ic_s,
        moc_cutoff=m.moc_cutoff,
        moc_beta_wdr=m.moc_beta_wdr,
        moc_offset_wdr=m.moc_offset_wdr,
        moc_maxrate_wdr=m.moc_maxrate_wdr,
        moc_minrate_wdr=m.moc_minrate_wdr,
        moc_beta_ic=m.moc_beta_ic,
        moc_offset_ic=m.moc_offset_ic,
        moc_maxrate_ic=m.moc_maxrate_ic,
        moc_minrate_ic=m.moc_minrate_ic,
        moc_weight_wdr=m.moc_weight_wdr,
        moc_weight_ic=m.moc_weight_ic,
        moc_len_integ=m.moc_len_integ,
        dur_pad_left=0.05,
        clip_left=true,
        dur_pad_right=0.02,
        clip_right=false,
        fs=m.fs,
    )
end

@with_kw struct TwoFormantTone <: CarneyLabUtils2.Params
    # Params interface
    name::String="two_formant_tone"
    gitcommit::String=DrWatson.gitdescribe(projectdir())
    gitcommit_specific::Bool=false
    pathhead::String="\\home\\daniel\\cl_sim\\efferent"
    exclude_identifier_keys::Vector{Symbol}=Symbol[]
    display_keys::Vector{Symbol}=[]

    # Stimulus interface
    f0::Float64=100.0
    dur::Float64=0.3
    dur_ramp::Float64=0.01
    fs::Float64=100e3
    a1::Float64=1.0
    f1::Float64=1500.0
    b1::Float64=400.0
    a2::Float64=2.0
    f2::Float64=3000.0
    b2::Float64=400.0
    l::Float64=70.0
end

function CarneyLabUtils2.synthesize!(x::TwoFormantTone)
    @unpack f0, dur, dur_ramp, fs, a1, f1, b1, a2, f2, b2, l = x
    x = sum(map(f -> cosine_ramp(scale_dbspl(pure_tone(f, 0.0, dur, fs), 0.0), dur_ramp, fs), f0:f0:min(10e3, fs/2)))
    x = a1 .* two_pole_resonator(x, f1, b1, fs) .+ a2 .* two_pole_resonator(x, f2, b2, fs)
    x = scale_dbspl(x, l)
    return x
end