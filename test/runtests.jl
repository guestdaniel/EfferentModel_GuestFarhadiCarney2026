using Test
using AuditorySignalUtils
using AuditoryNerveFiber
using AuditoryMidbrain
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
    "c2", "ihc", "expon", "sout1", "syn", "hsr", "lsr", "cn"]
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

    # If we're looking at the synapse/rate, we need to get every 10th sample (this is because 
    # the original code used a linear interpolation to upsample back to the stimulus 
    # sampling rate, but the new code is actually simulated at 100 kHz, producing large
    # disparities between sample points)
    if stage in ["syn", "hsr", "lsr"]
        sim = sim[1:10:end]
    end

    # Return (possibly subsetted) data
    if stage == "control"
        # Control onset is messed up a bit because it starts out non-zero, so simply 
        # zero-padding old control signal isn't viable and we only want to look at the 
        # relevant pieces
        sim = sim[2000:end]
    elseif stage in ["sout1", "sout2", "syn", "hsr", "lsr"]
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
    new = sim_gfc2023_dict(x, f; dur_pad_left=0.0, clip_left=false)[stage]

    orig = postprocess_simulations(orig, stage, "zbc2014", f)
    new = postprocess_simulations(new, stage, "gfc2023", f)

    return orig, new
end

function run_2014_vs_2023_pure_tone(cf::Vector{Float64}, stage::String, f=1000.0, l=50.0)
    # Create stimulus 
    x = pt(f, l)

    # Simulate original and new responses
    orig = map(_cf -> sim_orig_dict(x, _cf)[stage], cf)
    new = sim_gfc2023_dict(x, cf; dur_pad_left=0.0, clip_left=false)[stage]

    # Postprocess all responses
    orig = map(x -> postprocess_simulations(x[1], stage, "zbc2014", x[2]), zip(orig, cf))
    new = map(x -> postprocess_simulations(x[1], stage, "gfc2023", x[2]), zip(new, cf))

    return orig, new
end

"""
    run_full_vs_simple_pure_tone(stage, f, l)

Simulates responses for full vs simplified model funcs at a stage for a short pure tone stimulus
"""
function run_full_vs_simple_pure_tone(stage::String, f=1000.0, l=50.0)
    # Crestimulus 
    x = pt(f, l)

    # Simulate original and new responses
    full = sim_gfc2023_dict(x, f; dur_pad_left=0.0, clip_left=false)[stage]
    simp = sim_gfc2023_wrapper_dict(x, f; dur_pad_left=0.0, clip_left=false)[stage]

    full = postprocess_simulations(full, stage, "gfc2023", f)
    simp = postprocess_simulations(simp, stage, "gfc2023", f)

    return full, simp
end

# ==========================================================================================
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# ~~~~ Configure and define
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# ==========================================================================================
# Stages to test, peripheral
# Note: we currently omit sout2, which is difficult to match exactly due to 
# resampling issues and small numerical discrepancies. Both outputs are still analyzed 
# visually in other testing code
stages_peripheral = ["control", "c1", "c2", "ihc", "expon", "sout1", "syn", "hsr", "lsr"]

# Stages to test, subcortical
stages_subcortical = ["cn", "ic"]

# Define parameter sets at which we'll test SFIE implementations
params_sfie = [
    (1.0e-3, 2e-3, 1e-3, 1.0, 0.5),  # τ_e, τ_i, d, a, s
    (1.5e-3, 1e-3, 2e-3, 2.0, 1.2),
]

# Define function to set relative tolerance for comparisons at each stage
function get_rtol(stage)
    if stage in ["sout1", "sout2", "syn", "hsr", "lsr"]
        return 0.10
    elseif stage == "gain"
        return 0.01
    else
        return 0.001
    end
end

# ==========================================================================================
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# ~~~~ Check SFIE implementation
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# ==========================================================================================
@testset "SFIE implementation" begin
    # ======================================================================================
    # Check that get_alpha_norm produces identical values in C implementation
    # ======================================================================================
    @testset "Computing filter coefficients for τ = $τ" for τ in [1e-3, 2e-3, 4e-3]
        @test begin
            b_julia, a_julia = AuditoryMidbrain.get_α_normalized(τ, 100e3, 1.0)
            b_c = zeros(2); a_c = zeros(3);
            ccall(
                (:get_alpha_norm, "C:\\Users\\dguest2\\cl_code\\Helios\\src\\model\\libgfc2023.so"),
                Cvoid,
                (
                    Cdouble,
                    Cdouble,
                    Cdouble,
                    Ptr{Cdouble},
                    Ptr{Cdouble},
                ),
                τ, 100e3, 1.0, b_c, a_c,
            )
            (b_c ≈ b_julia) & (a_c ≈ a_julia)
        end
    end

    # ======================================================================================
    # Check that filter_alpha produces identical values in C implementation
    # ======================================================================================
    @test begin
        # Julia version
        b, a = AuditoryMidbrain.get_α_normalized(1e-3, 100e3, 1.0)
        x = sin.(2π .* 250.0 .* (0.0:(1/100e3):(0.1 - 1/100e3)))
        y_julia = filt(b, a, x)

        # C version
        b, a = AuditoryMidbrain.get_α_normalized(1e-3, 100e3, 1.0)
        x = sin.(2π .* 250.0 .* (0.0:(1/100e3):(0.1 - 1/100e3)))
        y_c = zeros(length(x))
        for i = 0:1:(length(x)-1)
            ccall(
                (:filter_alpha, "C:\\Users\\dguest2\\cl_code\\Helios\\src\\model\\libgfc2023.so"),
                Cvoid,
                (
                    Ptr{Cdouble},
                    Cint,
                    Cdouble,
                    Ptr{Cdouble},
                    Ptr{Cdouble},
                    Ptr{Cdouble},
                ),
                x, i, 100e3, b, a, y_c,
            )
        end

        # Compare
        y_c ≈ y_julia
    end
end

# ==========================================================================================
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# ~~~~ Check wrapper features
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# ==========================================================================================
# @testset "Wrapper features" begin
#     # ======================================================================================
#     # Check that right clipping doesnt affect results
#     # ======================================================================================
#     @test begin
#         unpadded = sim_gfc2023_dict(pt(), 1000.0; dur_pad_left=0.0, dur_pad_right=0.0, clip_right=false)["hsr"]
#         padded = sim_gfc2023_dict(pt(), 1000.0; dur_pad_left=0.0, dur_pad_right=0.1, clip_right=true)["hsr"]
#         all(unpadded .== padded)
#     end
# end

# ==========================================================================================
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# ~~~~ Check whether full and simple model functions yield same outputs
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# ==========================================================================================
# @testset "Full vs wrapper model" begin
#     @testset "1-kHz pure tone, 50 dB SPL, stage: $stage" for stage in ["ihc", "hsr", "lsr", "ic"]
#         full, simp = run_full_vs_simple_pure_tone(stage)
#         @test isapprox(full, simp; rtol=get_rtol(stage))
#     end
# end

# ==========================================================================================
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# ~~~~ Check whether new model outputs match 2014 model outputs (single channel)
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# ==========================================================================================
@testset "Regression vs 2014 --- single channel" begin
    # ======================================================================================
    # Check response to 1 kHz pure tone at every filter output
    # ======================================================================================
    @testset "1-kHz pure tone, 50 dB SPL, stage: $stage" for stage in stages_peripheral
        orig, new = run_2014_vs_2023_pure_tone(stage)
        @test isapprox(orig, new; rtol=get_rtol(stage))
    end
end

# ==========================================================================================
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# ~~~~ Check whether new model outputs match 2014 model outputs (multichannel)
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# ==========================================================================================
# @testset "Regression vs 2014 --- multichannel" begin
#     # ======================================================================================
#     # Check response to 1 kHz pure tone at 1 kHz and 2 kHz CFs
#     # ======================================================================================
#     @testset "1-kHz pure tone, 50 dB SPL, multichannel, stage: $stage" for stage in stages_peripheral
#         orig, new = run_2014_vs_2023_pure_tone([1000.0, 2000.0], stage)
#         pairs = zip(orig, new)
#         @test all(map(pair -> isapprox(pair[1], pair[2]; rtol=get_rtol(stage)), pairs))
#     end
# end

# ==========================================================================================
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# ~~~~ Check whether subcortial model outputs look reasonable and are behaving
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# ==========================================================================================
@testset "Regression vs 2004 --- single channel" begin
    # ======================================================================================
    # Check response to 1 kHz pure tone at cochlear nucleus
    # ======================================================================================
    @testset "CN at parameter values: $params" for params in params_sfie
        # Extract params
        τ_e, τ_i, d, a, s = params

        # Create stimulus
        x = pt(1000.0, 50.0)

        # Simulate response from C subcortical model
        out = sim_gfc2023_dict(
            x, 
            1000.0;
            cn_tau_e=τ_e,
            cn_tau_i=τ_i,
            cn_delay=d,
            cn_amp=a,
            cn_inh=s,
            dur_pad_left=0.0,
            clip_left=false,
        )
        new = out["cn"]

        # Simulate response from AuditoryMidbrain.jl for cochlear nucleus stage
        old = sim_sfie_nc2004(
            out["hsr"], 
            τ_e=τ_e,
            τ_i=τ_i,
            d_i=d,
            S=s,
            A=a,
        )

        # Compare
        @test old ≈ new
    end
    # ======================================================================================
    # Check response to 1 kHz pure tone at inferior colliculus
    # ======================================================================================
    @testset "IC at parameter values: $params" for params in params_sfie
        # Extract params
        τ_e, τ_i, d, a, s = params

        # Create stimulus
        x = pt(1000.0, 50.0)

        # Simulate response from C subcortical model
        out = sim_gfc2023_dict(
            x, 
            1000.0;
            ic_tau_e=τ_e,
            ic_tau_i=τ_i,
            ic_delay=d,
            ic_amp=a,
            ic_inh=s,
            dur_pad_left=0.0,
            clip_left=false,
        )
        new = out["ic"]

        # Simulate response from AuditoryMidbrain.jl for cochlear nucleus stage
        old = sim_sfie_nc2004(
            out["hsr"], 
            τ_e=0.5e-3,
            τ_i=2.0e-3,
            d_i=1.0e-3,
            S=0.6,
            A=1.5,
        )
        old = sim_sfie_nc2004(
            old,
            τ_e=τ_e,
            τ_i=τ_i,
            d_i=d,
            S=s,
            A=a,
        )

        # Compare
        @test old ≈ new
    end
end

@testset "Regression vs 2004 --- multichannel" begin
    # ======================================================================================
    # Check response to 1 kHz pure tone at inferior colliculus
    # ======================================================================================
    @testset "CN at parameter values: $params" for params in params_sfie
        # Extract params
        τ_e, τ_i, d, a, s = params

        # Create stimulus
        x = pt(1000.0, 50.0)

        # Simulate response from C subcortical model
        out = sim_gfc2023_dict(
            x, 
            [1000.0, 2000.0];
            ic_tau_e=τ_e,
            ic_tau_i=τ_i,
            ic_delay=d,
            ic_amp=a,
            ic_inh=s,
            dur_pad_left=0.0,
            clip_left=false,
        )
        new = out["ic"]

        # Simulate response from AuditoryMidbrain.jl for cochlear nucleus stage
        old = map(out["hsr"]) do x 
            cn = sim_sfie_nc2004(
                x, 
                τ_e=0.5e-3,
                τ_i=2.0e-3,
                d_i=1.0e-3,
                S=0.6,
                A=1.5,
            )
            ic = sim_sfie_nc2004(
                cn,
                τ_e=τ_e,
                τ_i=τ_i,
                d_i=d,
                S=s,
                A=a,
            )
            return ic
        end

        # Compare old to new
        pairs = zip(old, new)
        @test all(map(pair -> isapprox(pair[1], pair[2]), pairs))
    end
end

# ==========================================================================================
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# ~~~~ Check whether Julia and Mex wrappers provide same outputs
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# ==========================================================================================
# @testset "Julia vs Mex" begin
#     # First, we'll compare the Julia wrapper to the Mex wrapper by simulating responses to a
#     # pure tone with gain control disabled and verifying that responses at each output stage
#     # available in the Mex wrapper (ihc, hsr, lsr, ic, and gain) produce matched outputs for
#     # a single-CF response
#     @testset "1-kHz pure tone, 50 dB SPL, gain control disabled, single channel, stage: $stage" for stage in ["ihc", "hsr", "lsr", "ic", "gain"]
#         # Synthesize pure tone
#         x = scale_dbspl(cosine_ramp(pure_tone(1000.0, 0.0, 0.3, 100e3), 0.01, 100e3), 50.0)

#         # Run both models with gain control disabled
#         julia = sim_gfc2023_dict(
#             x, 
#             1000.0; 
#             dur_pad_left=0.0, 
#             dur_pad_right=0.0,
#             moc_weight_ic=0.0,
#             moc_weight_wdr=0.0,
#         )[stage]
#         matlab = sim_gfc2023_wrapper_dict_mex(
#             x, 
#             1000.0;
#             moc_weight_ic=0.0,
#             moc_weight_wdr=0.0,
#         )[stage]
#         @test isapprox(julia, matlab; rtol=get_rtol(stage))
#     end

#     # Next, we'll repeat the same simulations above except that we will turn gain control 
#     # on with very typical parameter values
#     @testset "1-kHz pure tone, 50 dB SPL, gain control enabled, single channel, stage: $stage" for stage in ["ihc", "hsr", "lsr", "ic", "gain"]
#         # Synthesize pure tone
#         x = scale_dbspl(cosine_ramp(pure_tone(1000.0, 0.0, 0.3, 100e3), 0.01, 100e3), 50.0)

#         # Run both models with gain control disabled
#         julia = sim_gfc2023_dict(
#             x, 
#             1000.0; 
#             dur_pad_left=0.0, 
#             dur_pad_right=0.0,
#             moc_weight_ic=1.0,
#             moc_weight_wdr=1.0,
#         )[stage]
#         matlab = sim_gfc2023_wrapper_dict_mex(
#             x, 
#             1000.0;
#             moc_weight_ic=1.0,
#             moc_weight_wdr=1.0,
#         )[stage]
#         @test isapprox(julia, matlab; rtol=get_rtol(stage))
#     end
# end

