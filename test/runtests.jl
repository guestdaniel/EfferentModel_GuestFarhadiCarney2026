using Test
using AuditorySignalUtils
using AuditoryNerveFiber
using Statistics
using DSP
using Helios

# ==========================================================================================
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# ~~~~ Write functions to aid testing protocol
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# ==========================================================================================
"""
    postprocess_simulations(sim, stage, model, cf=1000.0)

Postprocesses a simulation based on which model produced it to allow for testing

# Arguments
- `sim`: Vector containing some simulation result for a given stage
- `stage`: String indicating which stage this stage corresponds to, from ["control", "c1", \
    "c2", "ihc", "expon", "sout1", "syn", "cn"]
- `model`: Which model the response came from, from ["zbc2014", "gfc2023"]
- `cf`: Characteristic frequency (Hz)
"""
function postprocess_simulations(sim, stage, model, cf=1000.0)
    # If we're looking at control, c1, or c2 in an old model, we need to shift signal by 
    # delaypoint samples to accomodate the fact that we shifted delay from immediately after
    # IHC in original code to immediately after middle ear filter in new code
    if (model == "zbc2014") & (stage in ["control", "c1", "c2"])
        delay = ccall(
            (:delay_cat, "C:\\Users\\dguest2\\cl_code\\Helios\\src\\model\\libgfc2023.so"),
            Cdouble,
            (
                Cdouble,
            ),
            cf,
        )
        delaypoint = Int(ceil(delay/(1/100e3)))
        sim = shiftsignal(sim, delaypoint)
    end

    # If we're looking at power-law synapse stages, we need to adjust for the fact that
    # the new model operates at a single continuous sampling rate, while the old model
    # operates at a lower internal sampling rate for the power-law synapse stage 
    # To account for this, we downsample the new outputs by simply selecting every
    # 10th sample. We also need to account for the fact that the original sout1 and 
    # sout2 (and the whole powerlaw in general) included zeropadding on the edges of
    # the IHC respons. This is eliminated in the new code and produces edge effects numerical
    # the temporal edge of the simulations.
    if (model == "gfc2023") & (stage in ["sout1", "sout2"])
        sim = sim[1:10:end]
    end
    if (model == "zbc2014") & (stage in ["sout1", "sout2"])
        delaypoint = Int(floor(7500/(cf/1e3))/10)
        sim = sim[(1:(length(sim) - delaypoint*2)) .+ delaypoint]
    end

    # If we're looking at the synapse, we need to get every 10th sample (this is because 
    # the original code used a linear interpolation to upsample back to the stimulus 
    # sampling rate, but the new code is actually simulated at 100 kHz, producing large
    # disparities between sample points)
    if stage == "syn"
        sim = sim[1:10:end]
    end

    # Return (possibly subsetted) data
    if stage == "control"
        # Control onset is messed up a bit because it starts out non-zero, so simply 
        # zero-padding old control signal isn't viable and we only want to look at the 
        # relevant pieces
        sim = sim[2000:end]
    elseif stage in ["sout1", "sout2", "syn"]
        # For synapse stuff, we need to avoid the initial few samples because the lack of
        # a "delaypoint" system in the new model creates onset irregularities
        sim = sim[50:2000]
    end

    return sim
end

"""
    run_2014_vs_2023_pure_tone(stage, f, l)

Simulates responses for old and new model at a stage for a short pure tone stimulus
"""
function run_2014_vs_2023_pure_tone(stage::String, f=1000.0, l=50.0)
    # Create stimulus 
    x = pt(f, l)

    # Simulate original and new responses
    orig = sim_orig_dict(x, f)[stage]
    new = sim_gfc2023_dict(x, f)[stage]

    orig = postprocess_simulations(orig, stage, "zbc2014", f)
    new = postprocess_simulations(new, stage, "gfc2023", f)

    return orig, new
end

function run_2014_vs_2023_pure_tone(cf::Vector{Float64}, stage::String, f=1000.0, l=50.0)
    # Create stimulus 
    x = pt(f, l)

    # Simulate original and new responses
    orig = map(_cf -> sim_orig_dict(x, _cf)[stage], cf)
    new = sim_gfc2023_dict(x, cf)[stage]

    # Postprocess all responses
    orig = map(x -> postprocess_simulations(x[1], stage, "zbc2014", x[2]), zip(orig, cf))
    new = map(x -> postprocess_simulations(x[1], stage, "gfc2023", x[2]), zip(new, cf))

    return orig, new
end

# ==========================================================================================
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# ~~~~ Configure and define
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# ==========================================================================================
# Stages to test
# Note: we currently omit sout2, which is difficult to match exactly due to 
# resampling issues and small numerical discrepancies. Both outputs are still analyzed 
# visually in other testing code
stages = ["control", "c1", "c2", "ihc", "expon", "sout1", "syn"]

# Define function to set relative tolerance for comparisons at each stage
function get_rtol(stage)
    if stage in ["sout1", "sout2", "syn"]
        return 0.10
    else
        return 0.001
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
        orig, new = run_2014_vs_2023_pure_tone(stage)
        @test isapprox(orig, new; rtol=get_rtol(stage))
    end
end

@testset "Regression vs 2014 --- multichannel" begin
    # ======================================================================================
    # Check response to 1 kHz pure tone at 1 kHz and 2 kHz CFs
    # ======================================================================================
    @testset "1-kHz pure tone, 50 dB SPL, multichannel, stage: $stage" for stage in stages
        orig, new = run_2014_vs_2023_pure_tone([1000.0, 2000.0], stage)
        pairs = zip(orig, new)
        @test all(map(pair -> isapprox(pair[1], pair[2]; rtol=get_rtol(stage)), pairs))
    end
end