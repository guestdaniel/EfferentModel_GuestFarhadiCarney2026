# This script tests whether the binaural components of the model code work as
# intended, including method dispatch in Julia to handle binaural audio input,
# checking return types, comparing model responses in the two ears and model
# responses in comparison to the binaural model.

# First, let's verify that both methods are callable
@testset "Binaural callable" begin
    # Verify that basic method is callable with monaural input
    sim_gfc2023(zeros(5000), [1e3])

    # Verify that basic method is callable with binaural input
    sim_gfc2023([zeros(5000), zeros(5000)], [1e3])
end

# Next, let's verify that the correct return types are provided
@testset "Binaural return types" begin
    # Verify that monaural method returns Vector{Vector{Vector{Float64}}}
    sim_gfc2023(zeros(5000), [1e3]) isa Vector{Vector{Vector{Float64}}}

    # Verify that binaural method returns Vector{Vector{Vector{Vector{Float64}}}}
    sim_gfc2023([zeros(5000), zeros(5000)], [1e3]) isa Vector{Vector{Vector{Vector{Float64}}}}
end

# Verify that for a random noise input that the monaural and binaural model 
# outputs all match by default; disable efferent gain control and ffGN for
# simplicity's sake
@testset "Binaural match to monaural regression" begin
    # Configure stimulus
    stim = randn(5000)
    
    # Run each model
    resp_monaural = sim_gfc2023(stim, [1e3]; moc_weight=[0.0], fractional=false)
    resp_binaural = sim_gfc2023([stim, stim], [1e3]; moc_weight=[0.0], fractional=false)

    # Check each waveform
    for i in eachindex(resp_monaural)
        @test isapprox(resp_monaural[i], resp_binaural[i][1])  # check left
        @test isapprox(resp_monaural[i], resp_binaural[i][2])  # check right
    end
end

# Verify that if we put two different stimuli in the two ears, the responses 
# are actually different
@testset "Binaural inputs" begin
    x1 = zeros(5000)
    x2 = randn(5000)    
    resp = sim_gfc2023([x1, x2], [1e3])
    for r in resp
        @test !isapprox(r[1], r[2])
    end
end