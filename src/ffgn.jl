# ##################################################################################################
# Fractional Gaussian noise synthesis
# ##################################################################################################
# Local version of caplab_zbc2014.ffGn_native that returns bit-identical output (for the same RNG
# state) but runs much faster for short stimuli. The noise is synthesized at a low rate and then
# upsampled by a large factor (10,000 at 100 kHz), which ffGn_native does with DSP.resample. That
# (a) designs a new 370,001-tap anti-imaging filter on every call and (b) computes every upsampled
# sample even though only the first N are kept. Here, the resampling filter is built once per
# upsampling factor and reused, and the filter stops after the first N output samples.

# Cache of resampling filters keyed by upsampling factor. Each filter holds state while in use and
# is shared by all callers, so this code is not safe to call from multiple threads at once.
const RESAMPLE_FILTERS = Dict{Int, FIRFilter{DSP.Filters.FIRInterpolator{Float64}, Float64}}()

"""
    upsample_first(y, resamp, n_out)

Upsamples `y` by the integer factor `resamp` and returns the first `n_out` samples.

Equivalent to `DSP.resample(y, resamp)[1:n_out]` (bit-identical), but reuses a cached filter and
only computes the requested samples.

Not thread-safe: the cached filter is shared and modified during each call, so do not call this
(or `ffGn_fast`) from multiple threads at once. Separate processes (e.g., Distributed) are fine.
"""
function upsample_first(y::Vector{Float64}, resamp::Int, n_out::Int)
    sf = get!(() -> FIRFilter(resample_filter(resamp), resamp), RESAMPLE_FILTERS, resamp)
    reset!(sf)

    # Compensate for the filter delay, as DSP.resample does
    τ = timedelay(sf)
    setphase!(sf, τ)

    # Filtering is causal, so the first n_out outputs only need the first n_in inputs of the
    # zero-padded signal
    n_in = ceil(Int, n_out / resamp + τ) + 2
    x = zeros(n_in)
    copyto!(x, 1, y, 1, min(n_in, length(y)))
    out = filt(sf, x)
    length(out) >= n_out || error("Upsampled noise is shorter than expected")
    return resize!(out, n_out)
end

"""
    ffGn_fast(N, tdres, Hinput, noiseType, mu, sigma; safety=4)
    ffGn_fast(N, tdres, Hinput, noiseType, mu)

Synthesizes `N` samples of fractional Gaussian noise for the AN synapse model.

Drop-in replacement for `caplab_zbc2014.ffGn_native` with the same arguments and bit-identical
output for the same RNG state. See `ffGn_native` for a description of the arguments.

Not thread-safe (see `upsample_first`); use separate processes rather than threads to run
simulations in parallel.
"""
function ffGn_fast(
    N::Int64,
    tdres::Float64,
    Hinput::Float64,
    noiseType::Float64,
    mu::Float64,
    sigma::Float64;
    safety::Int64=4
)
    # Start by handling noiseType
    if noiseType == 0.0
        return zeros(N)
    end

    # If noiseType != 0, we're synthesizing fractional Gaussian noise
    # First, we downsample the number of points
    resamp = Int(ceil(1e-1 / tdres))
    nop = N
    N = Int(ceil(N / resamp) + 1)
    if N < 10
        N = 10
    end

    # Next, determine if fGn or fBm should be produced
    if Hinput < 1.0
        H = Hinput
        fBn = false
    else
        H = Hinput - 1
        fBn = true
    end

    # Calculate fGn
    if H == 0.5
        y = randn(N)
    else
        Nfft = Int(2 ^ ceil(log2(2*(N-1))))
        NfftHalf = Int(round(Nfft / 2))

        k = [0:(NfftHalf-1); NfftHalf:-1:1]
        Zmag = 0.5 * ( (k.+1) .^ (2*H) -2*k .^ (2*H) + abs.(k .- 1) .^ (2*H))

        Zmag = real.(fft(Zmag))
        Zmag = sqrt.(Zmag)

        Z = Zmag .* (randn(Nfft) + randn(Nfft) .* 1im)

        y = real.(ifft(Z)) * sqrt(Nfft)

        y = y[1:(N+safety)]
    end

    # Convert fGn to fBn if needed
    if fBn
        y = cumsum(y)
    end

    # Resample to match AN model, keeping only the samples we need
    y = upsample_first(y, resamp, nop)

    return sigma .* y
end

function ffGn_fast(
    N::Int64,
    tdres::Float64,
    Hinput::Float64,
    noiseType::Float64,
    mu::Float64;
)
    if mu == 100.0
        ffGn_fast(N, tdres, Hinput, noiseType, mu, 200.0)
    elseif mu == 5.0
        ffGn_fast(N, tdres, Hinput, noiseType, mu, 30.0)
    elseif mu == 0.1
        ffGn_fast(N, tdres, Hinput, noiseType, mu, 3.0)
    else
        error("The specified mu parameter is not recognized, must be 100, 5, or 0.1!")
    end
end
