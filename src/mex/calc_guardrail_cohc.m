function cohcs = calc_guardrail_cohc(cfs, species)
% cohcs = CALC_GUARDRAIL_COHC(cfs, species)
%
% Calculates a "guardrail" COHC value at specified CF values.
% 
% Based on a series of simulations reported in [[cite]], guardrail COHC
% values are values of COHC for a specific CF and species that prevent ΔL
% and ΔTH values due to efferent activity from exceeding what is plausible
% given MOC electrical stimulation experiments. These simualations were
% parametrized in terms of normalized cochlear distance according to:
%     Greenwood, D. D. (1990). A cochlear frequency-position function for
%     several species—29 years later. The Journal of the Acoustical Society
%     of America, 87(6), 2592–2605. https://doi.org/10.1121/1.399052
%
% Hence, below we convert CF to a cochlear distance based on species and
% then determine an appropriate COHC value based on a pre-determined linear
% mapping between cochlear distance and COHC.

% Convert from CF to cochlear distance based on species
if species == 1  % cat
	cd = log10(cfs ./ 456.0 + 0.8) ./ 2.1;
elseif species == 2  % human
	cd = log10(cfs ./ 165.4 + 0.88) ./ 2.1;
else
	error("Species settings other than 1 (cat) or 2 (Human Shera) are not supported.");
end

% Map from cochlear distance to guardrail COHC value using "magic numbers"
% derived in [[cite]].
	cohcs = exp(-1.976321 + 1.767505 .* cd);
end

