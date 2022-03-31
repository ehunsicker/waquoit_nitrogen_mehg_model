function [MeHgWaterConc] = HgMethylation(N_conc,TotalHgConc)
%HGMETHYLATION calculates resulting MeHg concentration in water from N
%   Calculate fraction of mercury methylated from N
methyl_percent = 0.5 * 0.0313 * exp(0.1115 * N_conc);
%   Calculate MeHg concentration
MeHgWaterConc = TotalHgConc * methyl_percent;
end

