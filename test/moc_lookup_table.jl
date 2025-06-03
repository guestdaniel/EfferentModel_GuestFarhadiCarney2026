# This script tests whether the MOC lookup table system appropriately approximates the 
# underlying nonlinearity. This is tested by driving the model with a frozen noise input
# and all stochasticity disabled and then comparing the `gain` stage output with and without
# use of the LUT.
@testset "MOC LUT" begin
    @test begin
        # Sample random-noise waveform as stimulus, scale to 80 dB SPL to get broad driving
        x = scale_dbspl(randn(100_000), 80.0)

        # Extract gain output from single-channel model run with and without LUT
        r1 = sim_gfc2023_dict(x, [1e3]; moc_use_lut=false)["gain"][1]
        r2 = sim_gfc2023_dict(x, [1e3]; moc_use_lut=true)["gain"][1]

        isapprox(r1, r2; rtol=0.001)
    end
end
