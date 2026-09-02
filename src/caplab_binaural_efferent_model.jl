module caplab_binaural_efferent_model

using DrWatson
using DSP
using AuditorySignalUtils
using FFTW
using Match
using Libdl
using Profile
using Dates
using ZilanyBruceCarney2014

const libgfc2026 = joinpath(@__DIR__, "..", "deps", "model", "libgfc2026." * Libdl.dlext)

include("c_interface.jl")
include("wrappers.jl")
include("wrappers_orig.jl")
include("utils.jl")
include("test_utils.jl")

end
