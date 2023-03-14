using Helios

function testfunc(n=10)
    for _ in 1:n
        sim_gfc2023(zeros(5000), 1000.0)
    end
end

@time testfunc()

# ----- Informal tracking -----
# commit 5545bbc - 0.60s for 5000 samples