using CairoMakie
using Helios
using Helios_plots
using AuditorySignalUtils
using Dates
using BenchmarkTools
using CarneyLabUtils2
using GLM
using DataFrames
using Crayons

### Configure paths and settings
path_data = "\\home\\daniel\\cl_sim\\efferent\\benchmark"
path_fig = "\\home\\daniel\\cl_fig\\efferent\\benchmark"
fs = 100e3

### Figure 1: Plot execution time for single channel as number of samples increases
# Determine file names
fn_data = joinpath(path_data, string(today()) * "_1.jld2")
fn_fig = joinpath(path_fig, string(today()) * "_1.png")

# Estimate execution time
durs = LogRange(0.1, 5.0, 20)
if isfile(fn_data) & ~override
    df = load(fn_data)["df"]
else
    times = map(durs) do dur
        x = scale_dbspl(cosine_ramp(pure_tone(1000.0, 0.0, dur, fs), 0.01, fs), 50.0)
        @elapsed sim_gfc2023_dict(x, 1000.0)
    end
    df = DataFrame(dur=durs, time=times)
    save(fn_data, Dict("df" => df))
end

# Plot
set_theme!(theme_carney)
fig = Figure(; resolution=(300, 300))
ax = Axis(
    fig[1, 1], 
    xscale=log10, 
    yscale=log10, 
    xticks=filter(x -> (x >= minimum(df.dur)) & (x <= maximum(df.dur)), 10.0 .^ (-10:1:10)),
    xminorticks=IntervalsBetween(9),
    yticks=filter(x -> (x >= minimum(df.time)) & (x <= maximum(df.time)), 10.0 .^ (-10:1:10)),
    yminorticks=IntervalsBetween(9),
)
scatter!(ax, durs, times)
ax.xlabel = "Input duration (s)"
ax.ylabel = "Execution time (s)"
m = lm(@formula(time ~ dur), df)
df_hat = DataFrame(dur=LogRange(0.1, 5.0, 2000), time=zeros(2000))
ŷ = predict(m, df_hat)
lines!(ax, df_hat.dur, ŷ; color=:black)
Label(fig[0, :], "Input duration vs execution time"; tellwidth=false)
fig
save(fn_fig, fig)