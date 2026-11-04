export sim_spikes_zbc2014

"""
    sim_spikes_zbc2014(synout::Vector{Float64}; fs=100e3, rng=default_rng())
    sim_spikes_zbc2014(synout::Vector{Vector{Float64}}; fs=100e3, rng=default_rng())
    sim_spikes_zbc2014(synout::Vector{Vector{Vector{Float64}}}; fs=100e3, rng=default_rng())

Simulates spike times from synapse output using the Zilany, Bruce, and Carney (2014) spike
generator.

Generates one spike train for each channel of `synout`, which must be the synapse output
*before* refractory effects, i.e., `syn_hsr` or `syn_lsr` from `sim_gfc2023_dict` (or the
corresponding elements of the output of `sim_gfc2023`).

The nesting of `synout` determines the nesting of the output: a single time series
(`Vector{Float64}`) gives a single spike train (`Vector{Float64}`); monaural output
(`Vector{Vector{Float64}}`, channel × time) gives one spike train per channel; and binaural
output (`Vector{Vector{Vector{Float64}}}`, ear × channel × time) gives one spike train per
ear and channel. Spike times are in seconds, relative to the first sample of `synout`.

Each call generates one spike train per channel. To simulate repeated presentations, call
the function repeatedly, for example:
- `[sim_spikes_zbc2014(syn) for _ in 1:n_rep]` for spike trains from the same synapse
  output (with `fractional=true`, all repetitions then share one realization of the
  fractional Gaussian noise), or
- `[sim_spikes_zbc2014(sim_gfc2023_dict(x, cf; fractional=true)["syn_hsr"]) for _ in
  1:n_rep]` for spike trains from a new simulation for each repetition.

The spike generator's refractory state starts from an approximate initial condition at the
first sample of `synout`. By default, `sim_gfc2023` and `sim_gfc2023_dict` clip the
zero-padded time before the stimulus (`clip_left=true`), so spike trains generated from their
output start in this initial state at stimulus onset. To let the refractory state settle
before the stimulus, consider simulating with `clip_left=false`, then removing spikes that
occur before `dur_pad_left` and subtracting `dur_pad_left` from the remaining spike times.

# Arguments
- `synout`: synapse output before refractory effects (sp/s); see above for shapes.
- `fs::Float64`: sampling rate of `synout` (Hz); default 100 kHz.
- `rng::AbstractRNG`: random-number generator used to draw the spike trains; default
  `Random.default_rng()`. Pass a seeded generator (e.g., `Xoshiro(1)`) for reproducible
  spike trains.

# Returns
- Spike times (s), nested to match `synout`.
"""
function sim_spikes_zbc2014(
    synout::Vector{Float64};
    fs::Float64=100e3,
    rng::AbstractRNG=default_rng(),
)
    tdres = 1 / fs
    n_max = spikegen_zbc2014_max_spikes(length(synout), tdres)
    # The C code needs uniform random numbers in (0, 1); rand draws from [0, 1)
    rand_unif = rand(rng, n_max + 2)
    rand_unif .= 1.0 .- rand_unif
    sptime = zeros(n_max)
    n_spikes = spikegen_zbc2014!(sptime, synout, tdres, rand_unif)
    return resize!(sptime, n_spikes)
end

function sim_spikes_zbc2014(
    synout::Vector{<:Vector};
    fs::Float64=100e3,
    rng::AbstractRNG=default_rng(),
)
    return [sim_spikes_zbc2014(s; fs=fs, rng=rng) for s in synout]
end
