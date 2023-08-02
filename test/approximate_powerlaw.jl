@testset "Exponential process approximation to PLA" begin
    @testset "Sampling rate = $fs Hz" for fs in [50e3, 100e3, 200e3]
        cfs = [500.0, 1000.0, 2000.0, 4000.0]
        @testset "CF: $cf" for cf in cfs
            orig, new = run_2023_vs_2023(
                pt(cf, 50.0, 0.2, fs), 
                [cf],
                Dict{Symbol, Any}(:powerlaw_mode => 1, :fs => fs),
                Dict{Symbol, Any}(:powerlaw_mode => 2, :fs => fs),
            )
            @testset "1-kHz pure tone, 50 dB SPL, stage: $stage" for stage in testparams["stages_peripheral"] 
                @test isapprox(orig[stage][1], new[stage][1]; rtol=0.02)
            end
        end
    end
end
