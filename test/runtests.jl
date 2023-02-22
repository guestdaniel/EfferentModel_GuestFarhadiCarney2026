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
        meout = zeros(length(x))
        modelout = zeros(length(x))

        # Apply ME filter
        model!(x, 1000.0, 1, 1/100e3, length(x), 1.0, 1.0, 1, meout, modelout)

        # Verify output is non-zero
        DSP.rms(meout) != 0
        DSP.rms(modelout) != 0
    end

end
