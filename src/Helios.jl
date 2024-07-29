module Helios

using CairoMakie
using DrWatson
using DSP
using AuditoryNerveFiber
using AuditorySignalUtils
using FFTW
using Match
using Libdl
#using MATLAB  # for interop to check against MATLAB/Mex model implementation
using Profile
using Dates

include("c_interface.jl")
include("wrappers.jl")
#include("wrappers_matlab.jl")
include("utils.jl")
include("adaptation.jl")  # temporary file for analyzing approximate power law implnt
include("profile.jl") 

include("test_utils.jl")

end # module Helios

