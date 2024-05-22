function [ihcout, hsrout, lsrout, icout, gain] = sim_efferent_model(x, cf, args)
% SIM_EFFERENT_MODEL(x, cf) simulates efferent-model response to row-vector
% stimulus x at CFs in vector cf. 
% 
% Returned values are matrices of size (n_chan, n_sample), where 
% `n_sample=length(x)` and `n_chan=length(cf)`. The various output matrices 
% are described below:
%   1) Inner hair cell "voltage", in a.u.
%   2) High-spontaneous-rate auditory-nerve instantaneous rate, in sp/s
%   3) Low-spontaneous-rate auditory-nerve instantaneous rate, in sp/s
%   4) Inferior-colliculus (IC) band-enhanced rate, in sp/s
%   5) Time-varying cochlear gain factor (in [0, 1], where 0==no gain, 1==max gain)
%
% SIM_EFFERENT_MODEL(x, cf) passes evaluates the efferent model on the
% input sound-pressure waveform at particular CFs. 
%
% SIM_EFFERENT_MODEL(x, cf, species=2) runs the efferent model for a 
% species value of 2 (which corresponds to human tuning based on data from 
% Shera). Other model parameters, such as sampling rate or IC model
% parameters, are adjusted the same way by specifying a key-value
% combinations (e.g., ic_tau_e=0.5e-3 would set the inhibitory IC delay to
% 0.5 ms, or moc_cutoff=2.0 would set the MOC lowpass cutoff to 2 Hz). See
% below for more details about available parameters and their default
% values (which are always used unless otherwise specified).
%
% To examine a changelog, please see `README.txt`
%
% Arguments:
% - x: Vector containing input sound-pressure waveform (Pa)
%
% - cf: Vector containing characteristic frequencies for each channel 
%   in the simulation (Hz). For multi-channel simulations, we operate 
%   under the assumption that CFs are in order from lowest to highest and
%   are equidistant on a log-frequency scale.
%
% - args.fs: Sampling rate of the simulation (Hz). Note that inputs must be 
%   sampled at this sampling rate. 
%
% - args.cohc: Vector containing values for the Cohc parameter in each
%   channel (in [0, 1], where 0==no contribution of OHCs, 1==maximum 
%   contribution of OHCs)
%
% - args.cihc: Vector containing values for the Cihc parameter in each
%   channel (in [0, 1], where 0==zero amplitude from C1-path IHCs, 1==full
%   amplitude from C1-path IHCs)
%
% - args.species: Which species to simulate in the basilar membrane/inner 
%   hair cell stage, 1==cat, 2==human (Shera), 3== (Moore and Glasberg)
%
% - args.powerlaw_mode: Whether to use true power-law adaptation
%   (powerlaw_mode == 1) or an approximate power-law implementation using a
%   set of 100 parallel exponential adaptation processes with time
%   constants fit computationally to match true power-law adaptation
%   (powerlaw_mode == 2). 
%
% - args.cn_tau_e: Excitatory time constant in CN stage (s)
%
% - args.cn_tau_i: Inhibitory time constant in CN stage (s)
%
% - args.cn_delay: Inhibitory delay time in CN stage (s)
%
% - args.cn_amp: Excitatory strength in CN stage
%
% - args.cn_inh: Inhibitory strength in CN stage 
%
% - args.ic_tau_e: Excitatory time constant in IC stage (s)
%
% - args.ic_tau_i: Inhibitory time constant in IC stage (s)
%
% - args.ic_delay: Inhibitory delay time in IC stage (s)
%
% - args.ic_amp: Excitatory strength in IC stage
%
% - args.ic_inh: Inhibitory strength in IC stage 
%
% - args.moc_cutoff: Cutoff of the lowpass filter used in the MOC stage
%   (Hz). The default value of 0.64 Hz yields a filter that matches that
%   used in the older single-channel efferent model (i.e., it produces a
%   "decay constant" of exp(-2pi * 0.64/100e3) ~= 1-3.9998e-5, which 
%   matches the constant used in the "old efferent" code, see Farhadi et 
%   al. 2023).
%
% - args.moc_beta_wdr: "beta" parameter in the MOC input-output
%   nonlinearity for the wide-dynamic-range MOC pathway (a.u.)
%
% - args.moc_offset_wdr: "offset" parameter in the MOC input-output
%   nonlinearity for the wide-dynamic-range MOC pathway (a.u.)
%
% - args.moc_minrate_wdr: Minimum possible gain factor in MOC input-output
%   nonlinearity for the wide-dynamic-range MOC pathway
%
% - args.moc_beta_ic: "beta" parameter in the MOC input-output
%   nonlinearity for the IC MOC pathway (a.u.)
%
% - args.moc_offset_ic: "offset" parameter in the MOC input-output
%   nonlinearity for the IC MOC pathway (a.u.)
%
% - args.moc_weight_wdr: Scalar value multiplied with lowpass-filtered 
%   wide-dynamic-range MOC pathway signal before signal is passed through 
%   MOC input-output nonlinearity
%
% - args.moc_weight_ic: Scalar value multiplied with lowpass-filtered IC
%   MOC pathway signal before signal is passed through MOC input-output 
%   nonlinearity
%
% - args.moc_width_wdr: "Width" of the wide-dynamic-range cross-channel
%   "spread" (octaves). For example, a value of one octave means that each
%   wide-dynamic-range MOC signal will "spread" to all channels that have
%   CFs that fall within a band centered on the CF with a width of one
%   octave (within +/- one-half octave)
%
% - args.noiseType: Integer value determining whether we use empty matrices
%   (noiseType == -1), matrices of "frozen" fractional Gaussian noise 
%   (noiseType == 0), or matrices of "fresh" fractional Gaussian noise based
%   on the current global RNG state (noiseType == 1) as inputs for the
%   noise governing the stochastic behavior of the power-law synapse in the
%   auditory-nerve model. 
%
% - args.dur_settle: How long to simulate responses to silence before
%   simulating a response to input time-pressure waveform (s). Default of 
%   0.01 s (10 ms). This duration of time gives the various dynamic stages
%   of the model (e.g., AN adaptation, efferent gain control) a chance to
%   "settle in" to a more steady-state response regime before simulating
%   the stimulus response. If this is set to too short an interval, you may
%   see some weird response features at simulation onset.
% 
% - args.moc_delay: Delay time (s) between MOC responses and changes to
%   cochlear gain. The return value called `gain` does not reflect this
%   delay (i.e., it shows the time-varying gain factor BEFORE the delay is
%   applied). A delay time of, for example, 5 ms means that the peripheral
%   filter stage uses the gain factor from 5 ms in the past when
%   calculating filter outputs. The default value is set to 25 ms, but this
%   is subject to change as various sources of data are analyzed to
%   determine an appropriate value.
%
% - args.display_info: Displays information about currently selected
%   parameter values and model version to the console before running
%   the model. Useful for debugging.

arguments
    x (1,:)
    cf (1,:) {mustBeGreaterThanOrEqual(cf, 125.0), mustBeLessThanOrEqual(cf, 40e3)}
    args.fs {mustBeGreaterThanOrEqual(args.fs,50e3), mustBeLessThanOrEqual(args.fs,200e3)} = 100e3
    args.cohc {mustBeGreaterThanOrEqual(args.cohc, 0.0), mustBeLessThanOrEqual(args.cohc, 1.0)} = ones(size(cf))
    args.cihc {mustBeGreaterThanOrEqual(args.cihc, 0.0), mustBeLessThanOrEqual(args.cihc, 1.0)} = ones(size(cf))
	args.species {mustBeMember(args.species, [1, 2, 3])} = 2 
	args.powerlaw_mode {mustBeMember(args.powerlaw_mode, [1, 2])} = 2
	args.cn_tau_e {mustBeGreaterThan(args.cn_tau_e, 0.0)} = 0.5e-3
    args.cn_tau_i {mustBeGreaterThan(args.cn_tau_i, 0.0)} =  2.0e-3
    args.cn_delay {mustBeGreaterThanOrEqual(args.cn_delay, 0.0)} = 1.0e-3
    args.cn_amp {mustBeGreaterThan(args.cn_amp, 0.0)} = 1.5
    args.cn_inh {mustBeGreaterThanOrEqual(args.cn_inh, 0.0)} = 0.6
    args.ic_tau_e {mustBeGreaterThan(args.ic_tau_e, 0.0)} = 1.0/(10.0 * 64.0);  % BMF == 64 Hz
    args.ic_tau_i {mustBeGreaterThan(args.ic_tau_i, 0.0)} =  1.0/(10.0 * 64.0)*1.5
    args.ic_delay {mustBeGreaterThanOrEqual(args.ic_delay, 0.0)} = 1.0/(10.0 * 64.0)*2.0
    args.ic_amp {mustBeGreaterThan(args.ic_amp, 0.0)} = 1.0
    args.ic_inh {mustBeGreaterThanOrEqual(args.ic_inh, 0.0)} = 0.9
    args.moc_cutoff {mustBeGreaterThanOrEqual(args.moc_cutoff, 0.0)} = 0.64
    args.moc_beta_wdr {mustBeGreaterThanOrEqual(args.moc_beta_wdr, 0.0)} = 0.015;
    args.moc_offset_wdr {mustBeGreaterThanOrEqual(args.moc_offset_wdr, 0.0)} = 10*4.0; 
	args.moc_minrate_wdr = 0.1;
	args.moc_maxrate_wdr = 1.0;
    args.moc_beta_ic {mustBeGreaterThanOrEqual(args.moc_beta_ic, 0.0)} = 0.015;
    args.moc_offset_ic {mustBeGreaterThanOrEqual(args.moc_offset_ic, 0.0)} = 10*4.0;
	args.moc_minrate_ic = 0.1;
	args.moc_maxrate_ic = 1.0;
    args.moc_weight_wdr {mustBeGreaterThanOrEqual(args.moc_weight_wdr, 0.0)} = 4.0;
    args.moc_weight_ic {mustBeGreaterThanOrEqual(args.moc_weight_ic, 0.0)} = 4.0;
    args.moc_width_wdr {mustBeGreaterThanOrEqual(args.moc_width_wdr, 0.0)} = 0.5
    args.noiseType {mustBeMember(args.noiseType, [-1, 0, 1])} = 1
	args.display_info {mustBeMember(args.display_info, [0, 1])} = 0
	args.dur_settle {mustBeGreaterThanOrEqual(args.dur_settle, 0.0)} = 0.2;
	args.clip_settle {mustBeMember(args.clip_settle, [0, 1])} = 1
	args.moc_delay {mustBeGreaterThanOrEqual(args.moc_delay, 0.0)} = 0.025;
end

% Determine number of samples corresponding to settle time
dur_orig = length(x)/args.fs;
len_settle = round(args.dur_settle * args.fs);

% Determine number of channels and samples
n_chan = length(cf);
n_sample = length(x) + len_settle;

% Synthesize fractional Gaussian noise
if args.noiseType == -1     
	% Use matrix of zeros
    ffGn_lsr = zeros(n_chan, n_sample);
    ffGn_hsr = zeros(n_chan, n_sample);
else                    
	% Synthesize fGn noise
    ffGn_lsr = zeros(n_chan, n_sample);
    ffGn_hsr = zeros(n_chan, n_sample);
    for ii=1:n_chan
        ffGn_lsr(ii, :) = ffGn(n_sample, 1/args.fs, 0.9, args.noiseType, 0.1, 3.0);
        ffGn_hsr(ii, :) = ffGn(n_sample, 1/args.fs, 0.9, args.noiseType, 100.0, 200.0);
    end
end

% Zero pad stimulus according to settle time
x = [zeros(1, len_settle), x];

% If we want to use display_info, print now
if args.display_info
	% CF-related info
	n_cf = length(cf);

	% Identify species
	switch args.species
		case 1
			species_string = "Cat";
		case 2
			species_string = "Human (Shera/Oxenham tuning)";
		case 3
			species_string = "Human (Glasberg/Moore tuning)";
	end

	% Identify noise type
	switch args.noiseType
		case -1
			noiseType_string = "No fGN (not recommended)";
		case 0
			noiseType_string = "Frozen fGn (i.e., fixed seed)";
		case 1
			noiseType_string = "Normal fGn";
	end

	% Identify power-law mode
	switch args.powerlaw_mode
		case 1
			powerlaw_mode_string = "True PLA (warning, slow!)";
		case 2
			powerlaw_mode_string = "Approximate PLA";
	end

	% Determine COHC/CIHC values
	cohc_min = round(min(args.cohc), 2);
	cohc_max = round(max(args.cohc), 2);
	cihc_min = round(min(args.cihc), 2);
	cihc_max = round(max(args.cihc), 2);

	% Determine git information that is available
	[fp, ~, ~] = fileparts(which("sim_efferent_model"));
	[r, s] = system(sprintf("git -C %s describe --tags --first-parent --abbrev=7 --long --dirty --always", fp));

	% Get information about version number
	fid = fopen(fullfile(fp, "model.c"));
	fgetl(fid);
	line = fgetl(fid);
	fclose(fid);
	parts = strsplit(line, " ");
	version_number = parts{4};

	% Print
	tic;
	fprintf("=========================================================\n");
	fprintf("Running Carney lab efferent model\n");
	fprintf("Version %s, last updated 5/22/2024\n", version_number);
	fprintf("Running " + string(n_cf) + " channels with...\n")
	if r == 0
		fprintf("	Build:                   %s", s);
	end
	fprintf("	Path to MATLAB function: %s\n", which("sim_efferent_model"));
	fprintf("	Path to Mex function:    %s\n", which("sim_efferent_model_mex"));
	fprintf("	Species:                 " + species_string + "\n")
	fprintf("	fGn type:                " + noiseType_string + "\n")
	fprintf("	PLA implementation:      " + powerlaw_mode_string + "\n")
	if length(cf) == 1
		fprintf("	CF:                      %0.2f kHz\n", cf/1000)
		fprintf("	COHC:                    %0.2f\n", cohc);
		fprintf("	COHC:                    %0.2f\n", cihc);
	else
		fprintf("	CFs:                     %0.2f to %0.2f kHz\n", cf(1)/1000, cf(end)/1000)
		fprintf("	COHC:                    " + string(cohc_min) + " to " + string(cohc_max) + "\n")
		fprintf("	CIHC:                    " + string(cihc_min) + " to " + string(cihc_max) + "\n")
	end
end

% Call Mex wrapper for efferent model
[ihcout, hsrout, lsrout, icout, gain] = sim_efferent_model_mex( ...
    x, ...
    ffGn_hsr, ...
    ffGn_lsr, ...
    cf, ...
    n_chan, ...
    1/args.fs, ...
    args.cohc, ...
    args.cihc, ...
    args.species, ...
    args.ic_tau_e, ...
    args.ic_tau_i, ...
    args.ic_delay, ...
    args.ic_amp, ...
    args.ic_inh, ...
    args.moc_cutoff, ...
    args.moc_beta_wdr, ...
    args.moc_offset_wdr, ...
    args.moc_beta_ic, ...
    args.moc_offset_ic, ...
    args.moc_weight_wdr, ...
    args.moc_weight_ic, ...
    args.moc_width_wdr, ...
	args.powerlaw_mode, ...
	args.cn_tau_e, ...
	args.cn_tau_i, ...
	args.cn_delay, ...
	args.cn_amp, ...
	args.cn_inh, ...
	args.moc_minrate_wdr, ...
	args.moc_maxrate_wdr, ...
	args.moc_minrate_ic, ...
	args.moc_maxrate_ic, ...
	args.dur_settle, ...
	args.moc_delay ...
);

% Optionally show tic/toc results
if args.display_info
	elapsed = toc;
	if elapsed > 1
		fprintf("	Runtime:                 %0.2f s (%0.2f s/chan)\n", elapsed, elapsed/n_chan);
	else
		fprintf("	Runtime:                 %0.2f ms (%0.2f s/chan)\n", elapsed*1000, elapsed*1000/n_chan);
	end
	fprintf("	Runtime / sim time:      %0.2f (lower is better, calc. per channel)\n", elapsed/dur_orig/n_chan);
end

% Optionally remove settle time at beginning of simulations
if args.clip_settle
	ihcout = ihcout(:, (len_settle+1):end);
	hsrout = hsrout(:, (len_settle+1):end);
	lsrout = lsrout(:, (len_settle+1):end);
	icout = icout(:, (len_settle+1):end);
	gain = gain(:, (len_settle+1):end);
end
end