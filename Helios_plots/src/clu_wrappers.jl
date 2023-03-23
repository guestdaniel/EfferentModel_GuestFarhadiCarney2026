export AuditorySubcortexGFCv1

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
    fractional::Bool=true

    # Subcortical parameters
    cn_τₑ::Float64=0.5e-3
    cn_τᵢ::Float64=2.0e-3
    cn_d::Float64=1.0e-3
    cn_a::Float64=1.5
    cn_s::Float64=0.6
    ic_τₑ::Float64=1.0e-3
    ic_τᵢ::Float64=2.0e-3
    ic_d::Float64=1.0e-3
    ic_a::Float64=4.0
    ic_s::Float64=0.9
    
    # Efferent parameters
    moc_cutoff::Float64=0.2
    moc_beta::Float64=0.01
    moc_offset::Float64=0.0
    moc_maxrate::Float64=1.0
    moc_minrate::Float64=0.0
    moc_weight_wdr::Float64=0.0
    moc_weight_ic::Float64=0.0
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
        moc_beta=m.moc_beta,
        moc_offset=m.moc_offset,
        moc_maxrate=m.moc_maxrate,
        moc_minrate=m.moc_minrate,
        moc_weight_wdr=m.moc_weight_wdr,
        moc_weight_ic=m.moc_weight_ic,
        dur_pad_left=0.05,
        clip_left=true,
        dur_pad_right=0.02,
        clip_right=true,
    )[m.stage]
end

function get_full_response(m, x::AbstractVector{Float64})
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
        moc_beta=m.moc_beta,
        moc_offset=m.moc_offset,
        moc_maxrate=m.moc_maxrate,
        moc_minrate=m.moc_minrate,
        moc_weight_wdr=m.moc_weight_wdr,
        moc_weight_ic=m.moc_weight_ic,
        dur_pad_left=0.05,
        clip_left=true,
        dur_pad_right=0.02,
        clip_right=true,
    )[m.stage]
end