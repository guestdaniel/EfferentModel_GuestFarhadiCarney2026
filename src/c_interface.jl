export model!

function model!(
    px::Vector{Float64},
    ffGn_hsr::Vector{Vector{Float64}},
    ffGn_lsr::Vector{Vector{Float64}},
    cf::Vector{Float64},
    n_chan::Int64,
    tdres::Float64,
    totalstim::Int64, 
    cohc::Float64,
    cihc::Float64,
    species::Int64,
    spont::Float64,
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
)
    ccall(
            (:model, "C:\\Users\\dguest2\\cl_code\\Helios\\src\\model\\libgfc2023.so"), 
            Cvoid, # return type
            (
                Ptr{Cdouble}, # px
                Ptr{Ptr{Cdouble}}, # ffGn_hsr
                Ptr{Ptr{Cdouble}}, # ffGn_lsr
                Ptr{Cdouble}, # cf
                Cint,         # nchan
                Cdouble,      # tdres
                Cint,         # totalstim
                Cdouble,      # cohc
                Cdouble,      # cihc
                Cint,         # species
                Cdouble,      # spont
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
        )
end

