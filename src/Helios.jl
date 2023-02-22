module Helios

using DSP
using FFTW
const libgfc2023 = "C:\\Users\\dguest2\\cl_code\\Helios\\src\\model\\libgfc2023.so"

export middle_ear!, model!

function model!(
    px::Vector{Float64},
    cf::Float64,
    nrep::Int64,
    tdres::Float64,
    totalstim::Int64, 
    cohc::Float64,
    cihc::Float64,
    species::Int64,
    meout::Vector{Float64},
    modelout::Vector{Float64},
)
    ccall(
            (:model, libgfc2023),  # function call
            Cvoid,                 # return type
            (
                Ptr{Cdouble}, # px
                Cdouble,      # cf
                Cint,         # nrep
                Cdouble,      # tdres
                Cint,         # totalstim
                Cdouble,      # cohc
                Cint,         # cihc
                Cint,         # species
                Ptr{Cdouble}, # meout
                Ptr{Cdouble}, # modelout
            ),
            px, cf, nrep, tdres, totalstim, cohc, cihc, species, meout, modelout  # input args
        )
end



function middle_ear!(
    px::Vector{Float64},
    tdres::Float64,
    totalstim::Int64, 
    species::Int64,
    meout::Vector{Float64},
)
    ccall(
            (:middle_ear, libgfc2023),  # function call
            Cvoid,                      # return type
            (
                Ptr{Cdouble}, # px
                Cdouble,      # tdres
                Cint,         # totalstim
                Cint,         # species
                Ptr{Cdouble}, # meout
            ),
            px, tdres,totalstim, species, meout,  # input args
        )
end


end # module Helios

