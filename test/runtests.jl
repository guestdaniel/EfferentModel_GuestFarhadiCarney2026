using Test
using AuditorySignalUtils
using ZilanyBruceCarney2014
using AuditoryMidbrain
using Statistics
using DSP
using EfferentModel_GuestFarhadiCarney2026

# Test implementation of SFIE in C
#include("sfie_implementation.jl")

# Test regression against 2014 model in single-channel simulations
include("regression_2014_single_channel.jl")

# Test that multichannel model results are identical to single-channel model results
include("singlechannel_vs_multichannel.jl")

# Test differences between true and approximate power-law adaptation
include("approximate_powerlaw.jl")

# Test implementation of normal PDF function in C
include("normal_pdf.jl")