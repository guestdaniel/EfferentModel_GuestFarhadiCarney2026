using Helios
using Profile
using ProfileView
using BenchmarkTools

# ##########################################################################################
# Case #1: MOC input-output nonlinearity overhead
# 
# The MOC input-output nonlinearity is:
#   C:    f(x) = (maxrate-minrate) * 1.0/(1.0 + pow(beta*(x-offset), 2.0))) + minrate
#   Math: f(x) = (maxrate-minrate) / (1 + (β*(x-offset))^2) + minrate
#
# This computation ends up being quite expensive because it has to be computed for every
# time step at the end of the time loop, and for every channel. If offset is set to a very 
# high value, the branch in `moc_nonlinearity` results in most time steps evaluating 
# instantly to 1, but if offset is set to a lower value or 0, then the nonlinearity must 
# be computed on every time step becasue x always > offset.
#
# Solution:
# The best strategy here would seem to be 
# 
# History:
#   - 6/2/2025 (prefix), 5 channel / 1e5 time samples / 5 reps 
#       moc_offset = 0.0     ->  7.236 s
#       moc_offset = 100000  ->  5.426 s
#
#   - 6/3/2025 (postfix), 5 channel / 1e5 time samples / 5 reps 
#       moc_offset = 0.0        ->  7.076 s
#       moc_offset = 100000     ->  6.076 s
#       moc_offset = 0.0 w/ lut ->  6.700 s
# ##########################################################################################
function test(; N=10_000, n_rep=5, moc_offset=0.0, moc_use_lut=false)
    x = scale_dbspl(randn(N), 80.0)
    for _ in 1:n_rep
        sim_gfc2023_dict(x, [1e3, 1e3, 1e3, 1e3, 1e3]; moc_offset=moc_offset, moc_use_lut=moc_use_lut)
    end
end

@btime test(; moc_offset=0.0, moc_use_lut=false)
@btime test(; moc_offset=0.0, moc_use_lut=true)
@btime test(; moc_offset=100000.0, moc_use_lut=false)
@btime test(; moc_offset=100000.0, moc_use_lut=true)

# The bulkiest three functions int 