module Helios

using DSP
using AuditoryNerveFiber
using AuditorySignalUtils
using FFTW
using MATLAB  # for interop to check against MATLAB/Mex model implementation

include("c_interface.jl")
include("wrappers.jl")
include("wrappers_matlab.jl")
include("utils.jl")

end # module Helios

