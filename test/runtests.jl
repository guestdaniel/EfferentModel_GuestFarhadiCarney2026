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
pt(f=1000.0, l=50.0, dur=0.2, fs=100e3) = scale_dbspl(pure_tone(f, 0.0, dur, fs), l)
stages = ["control", "c1", "c2", "ihc", "expon", "powerlaw", "sout1", "sout2", "syn"]

# ==========================================================================================
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# ~~~~ Check whether bindings are callable
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# ==========================================================================================
@testset "C bindings: callable" begin
    # ======================================================================================
    # model
    # ======================================================================================
    @test begin
        x = pt()
        cf = 1000.0
        meout_new = zeros(length(x))
        controlout_new = zeros(length(x))
        c1out_new = zeros(length(x))
        c1vihcout_new = zeros(length(x))
        c2out_new = zeros(length(x))
        c2vihcout_new = zeros(length(x))
        ihcout_new = zeros(length(x))
        synout_new = zeros(length(x))
        exponout_new = zeros(length(x))
        delaypoint = Int(floor(7500 / (cf / 1e3)))
        powerlawin_new = zeros(length(x) + delaypoint*3)
        fs = 100e3
        len_noise = Int(ceil((length(ihcout_new) + 2 * floor(7500 / (cf / 1e3))) * 1/fs * 10e3))
        ffGn = zeros(len_noise)
        sout1 = zeros(Int(ceil((length(x)+2*delaypoint) * 1/100e3 * 10e3)))
        sout2 = zeros(Int(ceil((length(x)+2*delaypoint) * 1/100e3 * 10e3)))

        # Run model
        model!(
            x, 
            ffGn,
            1000.0, 
            1/100e3, 
            length(x), 
            1.0, 
            1.0, 
            2, 
            100.0,
            0.0,
            1.0,
            meout_new, 
            controlout_new, 
            c1out_new, 
            c1vihcout_new, 
            c2out_new, 
            c2vihcout_new, 
            ihcout_new,
            synout_new,
            exponout_new,
            powerlawin_new,
            sout1,
            sout2,
        )

        true
    end
end

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

        # Calculate delaypoint if needed
        delay = ccall(
            (:delay_cat, "C:\\Users\\dguest2\\cl_code\\Helios\\src\\model\\libgfc2023.so"),
            Cdouble,
            (
                Cdouble,
            ),
            1000.0
        )
        delaypoint = Int(ceil(delay/(1/100e3)))

        # If we're looking at control, c1, or c2, we need to shift signal by delaypoint
        # samples to accomodate the fact that we shifted delay from immediately after IHC
        # in original code to immediately after middle ear filter in new code
        if stage in ["control", "c1", "c2"]
            orig = shiftsignal(orig, delaypoint)
        end

        # Verify fidelity using isapprox
        # TODO: determine appropriate tolerance values?
        if stage == "control"
            # Control onset is messed up a bit because it starts out non-zero, so simply 
            # zero-padding old control signal isn't viable
            @test orig[2000:end] ≈ new[2000:end]
        else
            @test orig ≈ new
        end
    end
end