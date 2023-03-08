export model!

function model!(
    px::Vector{Float64},
    ffGn::Vector{Vector{Float64}},
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
    exponout::Vector{Vector{Float64}},
    sout1::Vector{Vector{Float64}},
    sout2::Vector{Vector{Float64}},
    synout::Vector{Vector{Float64}},
)
    ccall(
            (:model, "C:\\Users\\dguest2\\cl_code\\Helios\\src\\model\\libgfc2023.so"), 
            Cvoid, # return type
            (
                Ptr{Cdouble}, # px
                Ptr{Ptr{Cdouble}}, # ffGn
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
                Ptr{Ptr{Cdouble}}, # exponOut
                Ptr{Ptr{Cdouble}}, # sout1
                Ptr{Ptr{Cdouble}}, # sout2
                Ptr{Ptr{Cdouble}}, # synout
            ),
            px,
            ffGn,
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
            exponout,
            sout1,
            sout2,
            synout,
        )
end

