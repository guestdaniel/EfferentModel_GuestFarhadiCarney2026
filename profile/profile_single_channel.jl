using Helios
using Profile
using ProfileView

# Profile overall model call
function test()
    sim_gfc2023_dict(zeros(500_000), 1e3)
end
Profile.clear()
@profile test()
ProfileView.view(; C=true)
Profile.print(; format=:flat, sortedby=:count, C=true)

# The bulkiest three functions int 