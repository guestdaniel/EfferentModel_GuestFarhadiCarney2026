using Helios

@time for rep in 1:10
    sim_anrate_gfc2023(zeros(5000), 1000.0)
end
@time sim_anrate_gfc2023(zeros(5000), repeat([1000.0], 10));

