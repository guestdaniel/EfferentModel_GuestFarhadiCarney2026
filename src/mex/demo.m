%% Example #1: Complex tone of 6-10F0, 50 dB SPL per component, sine phase
% Simulation may take quite some time to run --- if it takes too long,
% try reducing the duration or the number of channels below

% Construct stimulus
fs = 100e3;                                      % sample rate (Hz)
fs_down = 10e3;                                  % sample rate for plot (Hz)
dur = 0.5;                                       % duration (seconds)
t = 0.0:(1/fs):(dur - 1/fs);                     % sample times (s)
n_cf = 41;                                       % number of channels (#)
cf = exp(linspace(log(500.0), log(10e3), n_cf)); % CFs (Hz)
f0 = 200.0;                                      % F0 (Hz)
x = zeros(size(t));
for harm_no = 6:10
	x = x + 20e-6 * 10^(50.0/20.0) * sin(2*pi* harm_no*f0 * t)/0.7071;
end

% Simulate efferent-model response
% Note: If you want to change efferent-model parameters, this is the place
% to do it --- use the new key-value syntax (e.g., parametername=value) to
% specify values for different parameters. See documentation inside
% SIM_EFFERENT_MODEL for info about each parameter
[~, hsr, lsr, ic, gain] = sim_efferent_model(...
	x,...
	cf,...
	moc_weight_wdr=16.0,...
	moc_weight_ic=8.0,...
	moc_width_wdr=0.5...
);

% Resample responses down to lower sampling rate for plotting
t_resampled = 0.0:(1/fs_down):(dur - 1/fs_down);
hsr_resampled = resample(hsr, fs_down, fs, Dimension=2);
lsr_resampled = resample(lsr, fs_down, fs, Dimension=2);
ic_resampled = resample(ic, fs_down, fs, Dimension=2);
gain_resampled = resample(gain, fs_down, fs, Dimension=2);

% Plot as colorplot (HSR, LSR, IC, gain, in order)
figure;
tiledlayout(4, 1);
nexttile;
imagesc(t_resampled, 1:length(cf), hsr_resampled);
set(gca, 'ydir', 'normal');
caxis([0.0, 600.0]);
yticks(1:10:n_cf);
yticklabels(round(cf(1:10:n_cf)));
xlabel('Time (s)');
ylabel('CF (Hz)');
title('HSR');
nexttile;
imagesc(t_resampled, 1:length(cf), lsr_resampled);
set(gca, 'ydir', 'normal');
caxis([0.0, 200.0]);
yticks(1:10:n_cf);
yticklabels(round(cf(1:10:n_cf)));
xlabel('Time (s)');
ylabel('CF (Hz)');
title('LSR');
nexttile;
imagesc(t_resampled, 1:length(cf), ic_resampled);
set(gca, 'ydir', 'normal');
caxis([0.0, 100.0]);
yticks(1:10:n_cf);
yticklabels(round(cf(1:10:n_cf)));
xlabel('Time (s)');
ylabel('CF (Hz)');
title('IC (BE)');
nexttile;
imagesc(t_resampled, 1:length(cf), gain_resampled);
set(gca, 'ydir', 'normal');
caxis([0.0, 1.0]);
yticks(1:10:n_cf);
yticklabels(round(cf(1:10:n_cf)));
xlabel('Time (s)');
ylabel('CF (Hz)');
title('Gain');

%% Example #2: SAM noise MTF 
% Construct stimuli
fs = 100e3;                                      % sample rate (Hz)
dur = 0.5;                                       % duration (seconds)
t = 0.0:(1/fs):(dur - 1/fs);                     % sample times (s)
cf = 1000.0;                                     % CF (Hz)
fm_low = 2.0;                                    % lowest mod freq (Hz)
fm_high = 512.0;                                 % highest mod freq (Hz)
n_fm = 21;                                       % num mod freqs (#)
level = 10.0;                                    % spectrum level
depth = 0.0;                                     % modulation depth (dB)
m = 10^(depth/20);                               % modulation index
b_bp = fir1(4000, [100 8000]/(fs/2));            % filter coefs
fms = exp(linspace(log(fm_low), log(fm_high), n_fm));
stimuli = {};
for idx_fm = 1:n_fm
	% Construct carrier
	noise_spl = level+10*log10(fs/2);
	noise_rms = 10^(noise_spl/20)*20e-6;
    BBN_Pa = noise_rms*randn(1,dur*fs);
    BBN_Pa_band = conv(BBN_Pa,b_bp,'same');
    pre_mod_rms = rms(BBN_Pa_band);

	% Construct modulator
	modulator = m*sin(2*pi*fms(idx_fm)*t); 

	% Combine, scale, and store result
	temp = (1 + modulator) .* BBN_Pa_band;
	stimuli{idx_fm} = temp * pre_mod_rms / rms(temp);
end

% Get model responses
mu_with_eff = zeros(1, n_fm);
mu_wout_eff = zeros(1, n_fm);
for idx_fm = 1:n_fm
	% Call model w/ efferent system DISABLED
	[~, ~, ~, ic, ~] = sim_efferent_model(...
		stimuli{idx_fm},...
		cf,...
		moc_weight_wdr=0.0,...  % to disbale efferent system, set weights to 0
		moc_weight_ic=0.0...
	);
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

%% Example #3: Distribution of instantaneous spontaneous rates with ffGn settings
figure;

% Plot with "fresh" ffGn (noiseType == 1)
subplot(1, 2, 1);
hold on;
subtitle('Fresh ffGn (noiseType==1)');
for ii = 1:10
	[~, hsr, ~, ~, ~] = sim_efferent_model( ...
		zeros(1, 50000), ...
		[1000.0], ...
		noiseType=1 ...
	);
	t = 0.0:(1/100e3):(0.5-1/100e3);
	plot(t, hsr);
end
xlabel('Time (s)');
ylabel('Firing rate (sp/s)');
hold off;

% Plot with "frozen" ffGn (noiseType == 0)
subplot(1, 2, 2);
subtitle('Frozen ffGn (noiseType==0)');
hold on;
for ii = 1:10
	[~, hsr, ~, ~, ~] = sim_efferent_model( ...
		zeros(1, 50000), ...
		[1000.0], ...
		noiseType=0 ...
	);
	t = 0.0:(1/100e3):(0.5-1/100e3);
	plot(t, hsr);
end
hold off;
xlabel('Time (s)');
ylabel('Firing rate (sp/s)');

%% Example #4: Driven rates with fresh/frozen ffGn
figure;
fs = 100e3;                                      % sample rate (Hz)
dur = 0.5;                                       % duration (seconds)
t = 0.0:(1/fs):(dur - 1/fs);                     % sample times (s)
x = 20e-6 * 10^(50.0/20.0) * sin(2*pi * 1000.0 * t)*sqrt(2);

% Plot with "fresh" ffGn (noiseType == 1)
subplot(1, 2, 1);
hold on;
subtitle('Fresh ffGn (noiseType==1)');
for ii = 1:10
	[~, hsr, ~, ~, ~] = sim_efferent_model( ...
		x, ...
		[1000.0], ...
		noiseType=1 ...
	);
	plot(t, hsr);
end
xlabel('Time (s)');
ylabel('Firing rate (sp/s)');
hold off;

% Plot with "frozen" ffGn (noiseType == 0)
subplot(1, 2, 2);
subtitle('Frozen ffGn (noiseType==0)');
hold on;
for ii = 1:10
	[~, hsr, ~, ~, ~] = sim_efferent_model( ...
		x, ...
		[1000.0], ...
		noiseType=0 ...
	);
	plot(t, hsr);
end
hold off;
xlabel('Time (s)');
ylabel('Firing rate (sp/s)');