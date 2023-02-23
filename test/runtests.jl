using Test
using AuditorySignalUtils
using Statistics
using DSP
using Helios

pt(f=1000.0, l=50.0, dur=0.2, fs=100e3) = scale_dbspl(pure_tone(f, 0.0, dur, fs), l)

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
            1, 
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
