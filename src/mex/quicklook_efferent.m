function quicklook_efferent(stim, cf, fs, args)
	arguments
		stim
		cf
		fs=100e3
		args.moc_weight_ic=0.0
		args.moc_weight_wdr=0.0
		args.powerlaw_mode=2
	end

	% Pass stimulus through efferent model
	[ihc, hsr, lsr, ic, gain] = sim_efferent_model( ...
		[stim; zeros(5000, 1)], ...
		cf, ...
		noiseType=-1, ...
		moc_weight_ic=args.moc_weight_ic, ...
		moc_weight_wdr=args.moc_weight_wdr, ...
		powerlaw_mode=args.powerlaw_mode ...
	);

	% Plot everything
	figure;
	tiledlayout(5, 1, "TileSpacing", "compact", "Padding", "compact");
	responses = {ihc, hsr, lsr, ic, gain};
	labels = ["IHC", "HSR", "LSR", "IC", "Gain"];
	for idx_resp = 1:length(responses)
		resp = responses{idx_resp};
		t = 0.0:(1/fs):(length(resp)/fs - 1/fs);
		nexttile;
		plot(t, resp);
		title(labels(idx_resp));
		if idx_resp < length(responses)
			xticklabels([]);
		end
		if idx_resp == length(responses)
			ylim([0.0, 1.1]);
		else
			ylim([0.0, max(resp) * 1.2])
		end
	end
end

