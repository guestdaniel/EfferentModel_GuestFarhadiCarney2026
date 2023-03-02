using Test
using AuditorySignalUtils
using AuditoryNerveFiber
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
# Note: we currently omit synapse and sout2, which are difficult to match exactly due to 
# resampling issues and small numerical discrepancies. Both outputs are still analyzed 
# visually in other testing code
stages = ["control", "c1", "c2", "ihc", "expon", "sout1", "syn"]

# ==========================================================================================
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# ~~~~ Check whether new model outputs match 2014 model outputs
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# ==========================================================================================
@testset "Regression vs 2014" begin
    # ======================================================================================
    # Check response to 1 kHz pure tone at every filter output
    # ======================================================================================
    @testset "1-kHz pure tone, 50 dB SPL, stage: $stage" for stage in stages
        # Create stimulus 
        x = pt(1000.0, 50.0)

        # Simulate original and new responses
        orig = sim_orig_dict(x, 1000.0)[stage]
        new = sim_gfc2023_dict(x, 1000.0)[stage]

        # If we're looking at control, c1, or c2, we need to shift signal by delaypoint
        # samples to accomodate the fact that we shifted delay from immediately after IHC
        # in original code to immediately after middle ear filter in new code
        delay = ccall(
            (:delay_cat, "C:\\Users\\dguest2\\cl_code\\Helios\\src\\model\\libgfc2023.so"),
            Cdouble,
            (
                Cdouble,
            ),
            1000.0
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
            delaypoint = Int(floor(7500/(1000.0/1e3))/10)
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

        # Verify fidelity using isapprox
        # TODO: determine appropriate tolerance values?
        if stage == "control"
            # Control onset is messed up a bit because it starts out non-zero, so simply 
            # zero-padding old control signal isn't viable
            @test orig[2000:end] ≈ new[2000:end]
        elseif stage in ["sout1", "sout2", "syn"]
            # For synapse stuff, we need to avoid the initial few samples (which contain
            # an artifact due to the powerlaw "settling in") as well as use a more liberal 
            # relative tolerance, since I make no effort to precisely match the resampling 
            # strategies used in the model 
            @test orig[50:2000] ≈ new[50:2000] rtol=0.03
        else
            @test orig ≈ new
        end
    end
end