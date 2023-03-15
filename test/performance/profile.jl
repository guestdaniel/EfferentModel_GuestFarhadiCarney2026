using Helios
using Profile
using ProfileView

# Write function to be profiled
function profile_func(n=10)
    for _ in 1:n
        sim_gfc2023(zeros(5000), 1000.0)
    end
end

# Run once to trigger compilation
profile_func(1)

# Run profile and viz
Profile.clear()
@profile profile_func(100)
ProfileView.view(C=true)