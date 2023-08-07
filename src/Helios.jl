module Helios

using DSP
using AuditoryNerveFiber
using AuditorySignalUtils
using FFTW
using Match
#using MATLAB  # for interop to check against MATLAB/Mex model implementation

include("c_interface.jl")
include("wrappers.jl")
include("wrappers_ia.jl")
#include("wrappers_matlab.jl")
include("utils.jl")
include("adaptation.jl")  # temporary file for analyzing approximate power law implnt

include("test_utils.jl")

end # module Helios

