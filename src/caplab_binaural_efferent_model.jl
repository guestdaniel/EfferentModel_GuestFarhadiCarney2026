module caplab_binaural_efferent_model

using DSP
using Match
using Libdl

using caplab_signals
using caplab_zbc2014

const libgfc2026 = joinpath(@__DIR__, "..", "deps", "model", "libgfc2026." * Libdl.dlext)

include("c_interface.jl")
include("wrappers.jl")
include("wrappers_orig.jl")
include("utils.jl")
include("test_utils.jl")

__init__() = reset_model_library!()

end
