using Test
using AuditorySignalUtils
using Statistics
using DSP
using Helios

pt(f=1000.0, l=50.0, dur=0.2, fs=100e3) = scale_dbspl(pure_tone(f, 0.0, dur, fs), l)

# ==========================================================================================
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# ~~~~ Check whether bindings are callable and return non-zeros
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# ==========================================================================================
@testset "C bindings: callable" begin
    # ======================================================================================
    # middle_ear
    # ======================================================================================
    @test begin
        # Create inputs
        x = pt()
        meout = zeros(length(x))

        # Apply ME filter
        middle_ear!(x, 1/100e3, length(x), 1, meout)

        # Verify output is non-zero
        DSP.rms(meout) != 0
    end

    # ======================================================================================
    # model
    # ======================================================================================
    @test begin
        # Create inputs
        x = pt()
        meout_new = zeros(length(x))
        controlout_new = zeros(length(x))
        c1out_new = zeros(length(x))
        c1vihcout_new = zeros(length(x))
        c2out_new = zeros(length(x))
        c2vihcout_new = zeros(length(x))
        ihcout_new = zeros(length(x))

        # Run model
        model!(
            x, 
            1000.0, 
            1/100e3, 
            length(x), 
            1.0, 
            1.0, 
            2, 
            meout_new, 
            controlout_new, 
            c1out_new, 
            c1vihcout_new, 
            c2out_new, 
            c2vihcout_new, 
            ihcout_new
        )

        # Verify output is non-zero
        DSP.rms(ihcout_new) != 0
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
    @test begin
        # Run new model function
        x = scale_dbspl(pure_tone(1000.0, 0.0, 0.25, 100e3), 50.0)
        meout_new = zeros(length(x))
        controlout_new = zeros(length(x))
        c1out_new = zeros(length(x))
        c1vihcout_new = zeros(length(x))
        c2out_new = zeros(length(x))
        c2vihcout_new = zeros(length(x))
        ihcout_new = zeros(length(x))
        model!(
            x, 
            1000.0, 
            1/100e3, 
            length(x), 
            1.0, 
            1.0, 
            2, 
            meout_new, 
            controlout_new, 
            c1out_new, 
            c1vihcout_new, 
            c2out_new, 
            c2vihcout_new, 
            ihcout_new
        )

        # Call original model functions
        ihcout = zeros(length(x))
        c1out = zeros(length(x))
        c1vihcout = zeros(length(x))
        c2out = zeros(length(x))
        c2vihcout = zeros(length(x))
        controlout = zeros(length(x))
        ccall(
            (:IHCDEBUG, "C:\\Users\\dguest2\\cl_code\\Helios\\external\\julia\\libzbc2014debug.so"),
            Cvoid,                   # return type
            (                        # arg types
                Ptr{Cdouble},        # px
                Cdouble,             # cf
                Cint,                # nrep
                Cdouble,             # tdres
                Cint,                # totalstim
                Cdouble,             # cohc
                Cdouble,             # cihc
                Cint,                # species
                Ptr{Cdouble},        # ihcout
                Ptr{Cdouble},        # c1out
                Ptr{Cdouble},        # c1vihcout
                Ptr{Cdouble},        # c1out
                Ptr{Cdouble},        # c2vihcout
                Ptr{Cdouble},        # controlout
            ),
            x, 1000.0, 1, 1/100e3, length(x), 1.0, 1.0, 2, ihcout, c1out, c1vihcout, c2out, c2vihcout, controlout, # pass arguments
        )

        # Map through combos of outputs, check for identical responses
        matches = map(
            zip(
                [ihcout, controlout, c1out, c1vihcout, c2out, c2vihcout],
                [ihcout_new, controlout_new, c1out_new, c1vihcout_new, c2out_new, c2vihcout_new],
            )
        ) do (orig, new)
            all(orig .== new)
        end

        # Check that everything matches
        all(matches)

    end
end