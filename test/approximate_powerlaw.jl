@testset "Exponential process approximation to PLA" begin
    cfs = [500.0, 1000.0, 2000.0, 4000.0]
    orig, new = run_2023_vs_2023(
        pt(1000.0, 50.0), 
        cfs,
        Dict{Symbol, Any}(:powerlaw_mode => 1),
        Dict{Symbol, Any}(:powerlaw_mode => 2),
    )
    @testset "1-kHz pure tone, 50 dB SPL, stage: $stage" for stage in testparams["stages_peripheral"] 
        @testset "CF: $cf" for (idx, cf) in enumerate(cfs)
            @test isapprox(orig[stage][idx], new[stage][idx]; rtol=0.02)
        end
    end
end
