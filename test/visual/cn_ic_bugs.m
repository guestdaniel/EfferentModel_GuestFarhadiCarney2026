% Synthesize pure-tone stimulus
x = scale_dbspl(cosine_ramp(pure_tone(1000.0, 0.0, 0.2, 100e3), 0.01, 100e3), 50.0);

% Simulate HSR and LSR responses
hsr = sim_anrate_zbc2014(sim_ihc_zbc2014(x, 1000.0), 1000.0, fibertype=3, implnt=1);
lsr = sim_anrate_zbc2014(sim_ihc_zbc2014(x, 1000.0), 1000.0, fibertype=1, implnt=1);
[e_cn, i_cn] = filter_sfie_nc2004(hsr, 1e-3, 1e-3, 1e-3, 0.5, 1.0, 100e3);
cn_control = ((e_cn - i_cn) + abs(e_cn - i_cn))/2;
[e_ic, i_ic] = filter_sfie_nc2004(cn_control, 1e-3, 1e-3, 1e-3, 0.8, 1.0, 100e3);
ic_control = ((e_ic - i_ic) + abs(e_ic - i_ic))/2; 

% Visualize
subplot(4, 1, 1);
plot(hsr);
subplot(4, 1, 2);
plot(lsr);
subplot(4, 1, 3);
plot(cn_control);
subplot(4, 1, 4);
plot(ic_control);

% Simpler version
x = [1.0; zeros(9999, 1)];
[b, a] = get_alpha_norm(1e-3, 100e3, 1.0);
E = filter(b, a, x);
I = filter(b, a, x);
I = [zeros(fs*1e-3, 1); I];
I = I(1:length(E));
figure
hold on;
plot(I, 'b');
plot(E, 'r');
xlim([0.0, 1000.0]);
hold off;