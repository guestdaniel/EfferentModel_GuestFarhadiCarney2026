using Test
using AuditorySignalUtils
using AuditoryNerveFiber
using AuditoryMidbrain
using Statistics
using DSP
using Helios

# ==========================================================================================
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# ~~~~ Configure and define
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# ==========================================================================================
# Convenience function to synthesize pure tone
pt(f=1000.0, l=50.0, dur=0.2, fs=100e3) = scale_dbspl(pure_tone(f, 0.0, dur, fs), l)

# Stages to test
# Note: we currently omit sout2, which is difficult to match exactly due to 
# resampling issues and small numerical discrepancies. Both outputs are still analyzed 
# visually in other testing code
stages = ["control", "c1", "c2", "ihc", "expon", "sout1", "syn"]
function get_rtol(stage)
    if stage in ["sout1", "sout2", "syn"]
        return 0.10
    else
        return 0.001
    end
end

# Functions to run our simulations
function postprocess_simulations(orig, new, stage, cf=1000.0)
    # If we're looking at control, c1, or c2, we need to shift signal by delaypoint
    # samples to accomodate the fact that we shifted delay from immediately after IHC
    # in original code to immediately after middle ear filter in new code
    delay = ccall(
        (:delay_cat, "C:\\Users\\dguest2\\cl_code\\Helios\\src\\model\\libgfc2023.so"),
        Cdouble,
        (
            Cdouble,
        ),
        cf,
    )
    delaypoint = Int(ceil(delay/(1/100e3)))
    if stage in ["control", "c1", "c2"]
        orig = shiftsignal(orig, delaypoint)
    end

    # If we're looking at power-law synapse stages, we need to adjust for the fact that
    # the new model operates at a single continuous sampling rate, while the old model
    # operates at a lower internal sampling rate for the power-law synapse stage 
    # To account for this, we downsample the new outputs by simply selecting every
    # 10th sample. We also need to account for the fact that the original sout1 and 
    # sout2 (and the whole powerlaw in general) included zeropadding on the edges of
    # the IHC respons. This is eliminated in the new code.
    if stage in ["sout1", "sout2"]
        new = new[1:10:end]
        delaypoint = Int(floor(7500/(cf/1e3))/10)
        orig = orig[(1:length(new)) .+ delaypoint]
    end

    # If we're looking at the synapse, we need to get every 10th sample of both 
    # responses, because otherwise the test of equality will fail simply because the 
    # original code used a linear interpolation to upsample back to the stimulus 
    # sampling rate, but the new code is actually simulated at 100 kHz
    if stage in ["syn"]
        new = new[1:10:end]
        orig = orig[1:10:end]
    end

    # Return (possibly subset) data
    if stage == "control"
        # Control onset is messed up a bit because it starts out non-zero, so simply 
        # zero-padding old control signal isn't viable
        return orig[2000:end], new[2000:end]
    elseif stage in ["sout1", "sout2", "syn"]
        # For synapse stuff, we need to avoid the initial few samples (which contain
        # an artifact due to the powerlaw "settling in") as well as use a more liberal 
        # relative tolerance, since I make no effort to precisely match the resampling 
        # strategies used in the model 
        return orig[50:2000], new[50:2000]
    else
        return orig, new
    end
end

function test_run_simulations(stage::String)
    # Create stimulus 
    x = pt(1000.0, 50.0)

    # Simulate original and new responses
    orig = sim_orig_dict(x, 1000.0)[stage]
    new = sim_gfc2023_dict(x, 1000.0)[stage]

    postprocess_simulations(orig, new, stage)
end

function test_run_simulations(cf::Vector{Float64}, stage::String)
    # Create stimulus 
    x = pt(1000.0, 50.0)

    # Simulate original and new responses
    orig = map(_cf -> sim_orig_dict(x, _cf)[stage], cf)
    new = sim_gfc2023_dict(x, cf)[stage]

    # Postprocess all responses
    map(x -> postprocess_simulations(x[1], x[2], stage, x[3]), zip(orig, new, cf))
end

# ==========================================================================================
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# ~~~~ Check SFIE implementation
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# ==========================================================================================
@testset "SFIE implementation" begin
    # ======================================================================================
    # Check that get_alpha_norm produces identical values in C implementation
    # ======================================================================================
    @testset "Computing filter coefficients for τ = $τ" for τ in [1e-3, 2e-3, 4e-3]
        @test begin
            b_julia, a_julia = AuditoryMidbrain.get_α_normalized(τ, 100e3, 1.0)
            b_c = zeros(2); a_c = zeros(3);
            ccall(
                (:get_alpha_norm, "C:\\Users\\dguest2\\cl_code\\Helios\\src\\model\\libgfc2023.so"),
                Cvoid,
                (
                    Cdouble,
                    Cdouble,
                    Cdouble,
                    Ptr{Cdouble},
                    Ptr{Cdouble},
                ),
                τ, 100e3, 1.0, b_c, a_c,
            )
            (b_c ≈ b_julia) & (a_c ≈ a_julia)
        end
    end

    # ======================================================================================
    # Check that filter_alpha produces identical values in C implementation
    # ======================================================================================
    @test begin
        # Julia version
        b, a = AuditoryMidbrain.get_α_normalized(1e-3, 100e3, 1.0)
        x = sin.(2π .* 250.0 .* (0.0:(1/100e3):(0.1 - 1/100e3)))
        y_julia = filt(b, a, x)

        # C version
        b, a = AuditoryMidbrain.get_α_normalized(1e-3, 100e3, 1.0)
        x = sin.(2π .* 250.0 .* (0.0:(1/100e3):(0.1 - 1/100e3)))
        y_c = zeros(length(x))
        for i = 0:1:(length(x)-1)
            ccall(
                (:filter_alpha, "C:\\Users\\dguest2\\cl_code\\Helios\\src\\model\\libgfc2023.so"),
                Cvoid,
                (
                    Ptr{Cdouble},
                    Cint,
                    Cdouble,
                    Ptr{Cdouble},
                    Ptr{Cdouble},
                    Ptr{Cdouble},
                ),
                x, i, 100e3, b, a, y_c,
            )
        end

        # Compare
        y_c ≈ y_julia
    end
end

# ==========================================================================================
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# ~~~~ Check whether new model outputs match 2014 model outputs
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# ==========================================================================================
@testset "Regression vs 2014 --- single channel" begin
    # ======================================================================================
    # Check response to 1 kHz pure tone at every filter output
    # ======================================================================================
    @testset "1-kHz pure tone, 50 dB SPL, stage: $stage" for stage in stages
        orig, new = test_run_simulations(stage)
        @test isapprox(orig, new; rtol=get_rtol(stage))
    end
end

@testset "Regression vs 2014 --- multichannel" begin
    # ======================================================================================
    # Check response to 1 kHz pure tone at 1 kHz and 2 kHz CFs
    # ======================================================================================
    @testset "1-kHz pure tone, 50 dB SPL, multichannel, stage: $stage" for stage in stages
        pairs = test_run_simulations([1000.0, 2000.0], stage)
        @test all(map(pair -> isapprox(pair[1], pair[2]; rtol=get_rtol(stage)), pairs))
    end
end