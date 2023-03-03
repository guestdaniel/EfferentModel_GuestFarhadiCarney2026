export model!

function model!(
    px::Vector{Float64},
    ffGn::Vector{Float64},
    cf::Float64,
    tdres::Float64,
    totalstim::Int64, 
    cohc::Float64,
    cihc::Float64,
    species::Int64,
    spont::Float64,
    meout::Vector{Float64},
    controlout::Vector{Float64},
    c1out::Vector{Float64},
    c1vihcout::Vector{Float64},
    c2out::Vector{Float64},
    c2vihcout::Vector{Float64},
    ihcout::Vector{Float64},
    synout::Vector{Float64},
    exponout::Vector{Float64},
    sout1::Vector{Float64},
    sout2::Vector{Float64},
)
    ccall(
            (:model, "C:\\Users\\dguest2\\cl_code\\Helios\\src\\model\\libgfc2023.so"), 
            Cvoid, # return type
            (
                Ptr{Cdouble}, # px
                Ptr{Cdouble}, # ffGn
                Cdouble,      # cf
                Cdouble,      # tdres
                Cint,         # totalstim
                Cdouble,      # cohc
                Cdouble,      # cihc
                Cint,         # species
                Cdouble,      # spont
                Ptr{Cdouble}, # meout
                Ptr{Cdouble}, # control
                Ptr{Cdouble}, # c1 
                Ptr{Cdouble}, # c1vihc
                Ptr{Cdouble}, # c2 
                Ptr{Cdouble}, # c2vihc
                Ptr{Cdouble}, # ihcout
                Ptr{Cdouble}, # synout
                Ptr{Cdouble}, # exponOut
                Ptr{Cdouble}, # sout1
                Ptr{Cdouble}, # sout2
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

