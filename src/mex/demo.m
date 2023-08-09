%% Example #1: Complex tone of 6-10F0, 50 dB SPL per component, sine phase
% Here, we synthesize a five-component complex tone with frequencies at
% 6F0-10F0 for an F0 of 200 Hz at 50 dB SPL per component and then simulate 
% and visualize the efferent-model response to this stimulus.

% Set parameters
fs = 100e3;                                      % sample rate (Hz)
fs_down = 10e3;                                  % sample rate for plot (Hz)
dur = 0.5;                                       % duration (seconds)
f0 = 200.0;                                      % F0 (Hz)
t = 0.0:(1/fs):(dur - 1/fs);                     % sample times (s)  
n_cf = 41;                                       % number of channels (#)
cf = exp(linspace(log(500.0), log(3000.0), n_cf)); % CFs (Hz)

% Construct stimulus
x = zeros(size(t));
for harm_no = 6:10
	x = x + 20e-6 * 10^(50.0/20.0) * sin(2*pi* harm_no*f0 * t)/0.7071;
end

% Simulate efferent-model response
% Note: If you want to change efferent-model parameters, this is the place
% to do it --- use the new key-value syntax (e.g., parametername=value) to
% specify values for different parameters. See documentation inside
% SIM_EFFERENT_MODEL for info about each parameter. For example, here we
% specify the weights for WDR-driven and IC-driven gain control pathways, 
% the "bandwidth" of the WDR-driven gain control pathway, and which 
% powerlaw implementation to use.
[~, hsr, lsr, ic, gain] = sim_efferent_model(...
	x,...
	cf,...
	moc_weight_wdr=16.0,...
	moc_weight_ic=8.0,...
	moc_width_wdr=0.5,...
	powerlaw_mode=2 ...
);
responses = {hsr, lsr, ic, gain};  % store all output matrices in cell array

% Resample responses down to lower sampling rate for plotting
t_resampled = 0.0:(1/fs_down):(dur - 1/fs_down);
for ii = 1:4
	responses{ii} = resample(responses{ii}, fs_down, fs, Dimension=2);
end

% Define labels/limits for plot
labels = ["HSR", "LSR", "IC BE", "Gain"];
limits = {[0.0, 600.0], [0.0, 200.0], [0.0, 100.0], [0.0, 1.0]};

% Plot as colorplot
figure;
tiledlayout(4, 1);
for ii = 1:4
	nexttile;
	imagesc(t_resampled, 1:length(cf), responses{ii});
	set(gca, 'ydir', 'normal');
	caxis(limits{ii});
	yticks(1:10:n_cf);
	yticklabels(round(cf(1:10:n_cf)));
	xlabel('Time (s)');
	ylabel('CF (Hz)');
	title(labels(ii));
	colorbar;
end

%% Example #2: SAM noise MTF
% Here, we construct a series of sinusoidally amplitude-modulated (SAM) 
% noises with various modulation frequencies and measure the mean IC rate
% in response to each stimulus with and without efferent gain control. 
%
% The BE cell should have a range of modulation frequencies over which it
% exhibits an increased firing rate. At the sound level used (30 dB SPL 
% spectrum level), however, there is little modulation coding in the BE
% cell without efferent gain control. Enabling efferent gain control
% improves this coding and enhanced the BE MTF shape.

% Set parameters
fs = 100e3;                                      % sample rate (Hz)
dur = 0.5;                                       % duration (seconds)
t = 0.0:(1/fs):(dur - 1/fs);                     % sample times (s)
cf = 1000.0;                                     % CF (Hz)
fm_low = 2.0;                                    % lowest mod freq (Hz)
fm_high = 512.0;                                 % highest mod freq (Hz)
n_fm = 21;                                       % num mod freqs (#)
level = 30.0;                                    % spectrum level
depth = 0.0;                                     % modulation depth (dB)
m = 10^(depth/20);                               % modulation index
b_bp = fir1(4000, [100 8000]/(fs/2));            % filter coefs
fms = exp(linspace(log(fm_low), log(fm_high), n_fm));

% Construct stimuli
stimuli = {};
for idx_fm = 1:n_fm
	% Construct carrier
	noise_spl = level+10*log10(fs/2);
	noise_rms = 10^(noise_spl/20)*20e-6;
    BBN_Pa = noise_rms*randn(1,dur*fs);
    BBN_Pa_band = conv(BBN_Pa,b_bp,'same');  % bandpass filter carrier from 0.1-8 kHz
    pre_mod_rms = rms(BBN_Pa_band);

	% Construct modulator
	modulator = m*sin(2*pi*fms(idx_fm)*t); 

	% Combine, scale, and store result
	temp = (1 + modulator) .* BBN_Pa_band;
	stimuli{idx_fm} = temp * pre_mod_rms / rms(temp);
end

% Get model responses
mu_with_eff = zeros(1, n_fm);  % output IC means without efferent gain control
mu_wout_eff = zeros(1, n_fm);  % output IC means WITH efferent gain control
for idx_fm = 1:n_fm
	% Call model w/ efferent system DISABLED
	[~, ~, ~, ic, ~] = sim_efferent_model(...
		stimuli{idx_fm},...
		cf,...
		moc_weight_wdr=0.0,...  % to disbale efferent system, set WDR weight to 0
		moc_weight_ic=0.0...    % ... and also set IC weight to 0
	);
	% Average and store IC rate
	mu_wout_eff(idx_fm) = mean(ic);  % note: this includes onset!

	% Call model w/ efferent system enabled
	[~, ~, ~, ic, ~] = sim_efferent_model(...
		stimuli{idx_fm},...
		cf,...
		moc_weight_wdr=2.0,...
		moc_weight_ic=8.0...
	);
	
	% Average and store IC rate
	mu_with_eff(idx_fm) = mean(ic);  % note: this includes onset!
end

% Plot
figure;
plot(fms, mu_wout_eff, 'b'); hold on;
plot(fms, mu_with_eff, 'r'); hold off;
set(gca, 'xscale', 'log');
xlabel('Modulation frequency (Hz)');
ylabel('Firing rate (sp/s)');
legend(["Without efferent", "With efferent"])
ylim([0.0, 50.0]);

%% Example #3: Distributions of spontaneous rates
% In the 2009-and-beyond models, some amount of fractional Gaussian noise
% was added to the synapse output inside the power-law adaptation loop. 
% This noise has several effects, one of which is to cause the spontaneous
% rate of the model to "drift" over time with a long-range temporal
% dependence. 
%
% In some cases, this behavior is undesireable, so two options
% are provided for control. First, by passing `noiseType=0`, the "fresh"
% noise is replaced with a frozen sample of noise that is reused on
% subsequent calls to the model --- this makes the model completely
% deterministic. Second, by passing `noiseType=-1`, the noise is eliminated
% entirely from the model. This can have unintended side effects, so this
% is not recommended unless you know what you are doing (but is shown here
% anyway). 
% 
% Below, example spontaneous responses are shown for each noiseType
% setting. For clarity, when the "frozen" responses are plotted, a small 
% increment is added to each waveform to ensure they are not plotted 
% directly on top of each other. 

figure;
noiseTypes = [-1, 0, 1];
labels = [
	"No fGn (noiseType == -1)", ...
	"Frozen fGn (noiseType = 0)", ...
	"Fresh fGn (noiseType = 1)" ...
];

% Loop through possibilities, generate responses and plot 
for ii = 1:3
	% Set up plot and titles
	subplot(1, 3, ii);
	subtitle(labels(ii));
	xlabel('Time (s)');
	ylabel('Firing rate (sp/s)');
	ylim([0.0, 750.0]);

	% Run simulations 10 times
	hold on;
	for jj = 1:10
		% Run model
		[~, hsr, ~, ~, ~] = sim_efferent_model( ...
			zeros(1, 50000), ...
			[1000.0], ...
			noiseType=noiseTypes(ii) ...
		);

		% Handle offsetting
		if ii == 1 || ii == 2
			offset = jj;
		else
			offset = 0.0;
		end

		% Plot result
		t = 0.0:(1/100e3):(0.5-1/100e3);
		plot(t, hsr + offset);
	end
	hold off;
end

%% Example #4: Comparison of true vs approximate power law for many pure tones
% In the efferent model, one of the "hungriest" parts of the calculations
% is implementing power-law adaptation, due to its dependence on all
% previous time samples. One way to get around this is to replace "true"
% power-law adaptation with an approximation composed of many parallel
% exponential adaptation processes with time constants varying over a wide
% range (to capture the behavior of the power law at different time
% scales). These exponential adaptation processes can be efficiently
% implemented as IIR filters that depend only on the previous sample,
% rather than all previous samples (as for the "true" power-law
% adaptation). 
% 
% This demo simulates responses to pure tones under the "true"
% power-law adaptation or under its approximate implementation and compares
% them side-by-side.

% Set parameters
fs = 100e3;                                      % sample rate (Hz)
dur = 0.2;                                       % duration (s)
dur_post = 0.1;                                  % duration of post-stimulus simulation time (s)
t = 0.0:(1/fs):(dur+dur_post - 1/fs);            % sample times (s)
freqs = [500.0, 1000.0, 2000.0, 4000.0, 8000.0]; % test frequencies (Hz)
time_windows = {[0.0, 0.001], [0.05, 0.06], [0.1, 0.3]};  % time windows for analysis (s)

% Pre-allocate storage
true = zeros(length(t), length(freqs));
approx = zeros(length(t), length(freqs));

% Loop over stimuli and do calculations (in parallel using parfor)
parfor idx_freq = 1:length(freqs)
	% Synthesize stimulus for this frequency (50 dB SPL pure tone)
	stim = 20e-6 * 10^(50.0/20.0) * ...
		sin(2*pi * freqs(idx_freq) * (0.0:(1/fs):(dur - 1/fs))) * sqrt(2);
	stim = [stim zeros(1, round(dur_post*fs))];

	% Run efferent model with true and approximate power-law adaptation
	[~, true(:, idx_freq), ~, ~, ~] = ...
		sim_efferent_model(stim, freqs(idx_freq), powerlaw_mode=1, noiseType=0);
	[~, approx(:, idx_freq), ~, ~, ~] = ...
		sim_efferent_model(stim, freqs(idx_freq), powerlaw_mode=2, noiseType=0);	
end

% Create figure
figure;

% Loop through frequencies
for idx_freq = 1:length(freqs)
	% Loop through different time scales
	for idx_tw = 1:length(time_windows)
		% Create subplot
		subplot( ...
			length(freqs), ...
			length(time_windows), ...
			length(time_windows)*(idx_freq-1) + idx_tw ...
		);

		% Plot true and approximate results on top of each other
		plot(t, true(:, idx_freq)); hold on;
		plot(t, approx(:, idx_freq)); hold off;

		% Add legend to upper left
		if idx_freq == 1 && idx_tw == 1
			legend(["True power-law adaptation", "Approximate power-law adaptation"]);
		end

		% Set limits, titles, and labels
		xlim(time_windows{idx_tw});
		xlabel("Time (s)");
		ylabel("Firing rate (sp/s)");
		ylim([0.0, 1000.0]);
		title(sprintf( ...
			"Freq = %4.0f Hz, time window = [%4.3f, %4.3f] s", ...
			freqs(idx_freq), ...
			time_windows{idx_tw}(1), ...
			time_windows{idx_tw}(2) ...
		));
	end
end

%% Example #5: Performance benefits of approximate power-law implementation
% See above for more details about power-law adaptation approximation. 
%
% Here, we test the approximation and compare its performance to the true
% power-law adaptation implementation. Feel free to adjust the durations
% below to determine the performance gains you can expect for different
% stimuli, although beware that durations beyond ~1s start to require
% prohibitively long compute times for the true power law (so it may take 
% quite some time to generate). The first graph simply depicts compute
% time, while the other graph shows the "performance gain" as a time ratio 
% from switching to the approximate power-law adaptatation.
%
% Note that time estimates for shorter stimuli can be somewhat unreliable
% --- if the graphs look nonmonotonic or have outliers, do not trust those
% points, since they can be unduly influenced by brief delays or
% interruptions in compute due to system processes or other factors. Also
% since we estimate compute time multiple times to eliminate some of this
% randomness, the plot can take quite some time to generate.

% Set parameters
durs = exp(linspace(log(5e-2), log(5e-1), 20));  % stimulus durations (s)
n_rep = 5;                                       % how many repeats to do

% Pre-allocate storage
durs_true = zeros(length(durs), n_rep);
durs_approx = zeros(length(durs), n_rep);

% Time compute time for each stimulus duration
for ii = 1:length(durs)
	% Run model with real power-law adaptation
	for jj = 1:n_rep
		tic;
		sim_efferent_model(zeros(1, round(durs(ii)*100e3)), [1000.0], powerlaw_mode=1);
		durs_true(ii, jj) = toc;
	end

	% Run model with approximate power-law adaptation
	for jj = 1:n_rep
		tic;
		sim_efferent_model(zeros(1, round(durs(ii)*100e3)), [1000.0], powerlaw_mode=2);
		durs_approx(ii, jj) = toc;
	end
end

% Plot compute time
figure;
subplot(1, 2, 1);
errorbar(durs, mean(durs_true, 2), 1.96 * std(durs_true, 0, 2)/sqrt(n_rep), 'k'); hold on;
grid on;
errorbar(durs, mean(durs_approx, 2), 1.96 * std(durs_true, 0, 2)/sqrt(n_rep), 'r');
plot(10 .^ (-3.0:0.01:3.0), 10 .^ (-3.0:0.01:3.0), 'k'); hold off;  % plot unity line
set(gca, "xscale", "log");
set(gca, "yscale", "log");
ylim([1e-2, 1e2]);
xlim([1e-2, 1e1]);
legend(["True power law adaptation", "Approximate power law adaptation"]);
xlabel("Stimulus duration (s)");
ylabel("Compute time (s)");

% Plot speed up
subplot(1, 2, 2);
plot(durs,  mean(durs_true, 2) ./  mean(durs_approx, 2), 'k'); grid on;
set(gca, "xscale", "log");
set(gca, "yscale", "log");
ylim([1e0, 1e3]);
xlim([2e-2, 2e0]);
xlabel("Stimulus duration (s)");
ylabel("Performance gain (true/approx compute time)");
