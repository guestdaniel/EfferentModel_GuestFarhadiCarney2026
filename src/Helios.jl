module Helios

using DSP
using FFTW
const libgfc2023 = "C:\\Users\\dguest2\\cl_code\\Helios\\src\\model\\libgfc2023.so"

export middle_ear!

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

