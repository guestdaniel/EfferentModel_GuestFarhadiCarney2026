export model!

function model!(
    px::Vector{Float64},
    ffGn::Vector{Float64},
    cf::Vector{Float64},
    tdres::Float64,
    totalstim::Int64, 
    cohc::Float64,
    cihc::Float64,
    species::Int64,
    spont::Float64,
    meout::Vector{Float64},
    controlout::Vector{Vector{Float64}},
    c1out::Vector{Vector{Float64}},
    c1vihcout::Vector{Vector{Float64}},
    c2out::Vector{Vector{Float64}},
    c2vihcout::Vector{Vector{Float64}},
    ihcout::Vector{Vector{Float64}},
    synout::Vector{Vector{Float64}},
    exponout::Vector{Vector{Float64}},
    sout1::Vector{Vector{Float64}},
    sout2::Vector{Vector{Float64}},
)
    ccall(
            (:model, "C:\\Users\\dguest2\\cl_code\\Helios\\src\\model\\libgfc2023.so"), 
            Cvoid, # return type
            (
                Ptr{Cdouble}, # px
                Ptr{Cdouble}, # ffGn
                Ptr{Cdouble}, # cf
                Cdouble,      # tdres
                Cint,         # totalstim
                Cdouble,      # cohc
                Cdouble,      # cihc
                Cint,         # species
                Cdouble,      # spont
                Ptr{Cdouble}, # meout
                Ptr{Ptr{Cdouble}}, # control
                Ptr{Ptr{Cdouble}}, # c1 
                Ptr{Ptr{Cdouble}}, # c1vihc
                Ptr{Ptr{Cdouble}}, # c2 
                Ptr{Ptr{Cdouble}}, # c2vihc
                Ptr{Ptr{Cdouble}}, # ihcout
                Ptr{Ptr{Cdouble}}, # synout
                Ptr{Ptr{Cdouble}}, # exponOut
                Ptr{Ptr{Cdouble}}, # sout1
                Ptr{Ptr{Cdouble}}, # sout2
            ),
            px,
            ffGn,
            cf, 
            tdres, 
            totalstim, 
            cohc, 
            cihc, 
            species, 
            spont,
            meout, 
            controlout, 
            c1out, 
            c1vihcout,
            c2out, 
            c2vihcout,
            ihcout,
            synout,
            exponout,
            sout1,
            sout2,
        )
end

