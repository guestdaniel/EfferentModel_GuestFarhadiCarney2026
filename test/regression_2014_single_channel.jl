@testset "Regression vs 2014 --- single channel" begin
    # Loop over different frequencies
    @testset "CF/tone freq = $cf Hz" for cf in [1000.0, 2000.0, 4000.0]
        # Simulate responses (with zero-padding to emulate "delaypoint" system in new model)
        orig, new = run_2014_vs_2023(
            pt(cf, 50.0), 
            [cf],
            Dict{Symbol, Any}(),  # leave old param values at default values
            Dict{Symbol, Any}(:dur_pad_left => 7500/(cf/1e3)/100e3),
        )

        # Loop through each stage and verify match
        @testset "stage: $stage" for stage in testparams["stages_peripheral"] 
            @test isapprox(orig[stage], new[stage]; rtol=rtol_2014_vs_2023(stage)) broken=((cf == 4000.0) & (stage in ["syn", "lsr"]))
        end
    end
end
