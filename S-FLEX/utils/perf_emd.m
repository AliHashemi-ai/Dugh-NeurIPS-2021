function [em mean_co mean_co_dist in_co] = perf_emd(Weight, X_true, D, ind_true) 

amp_est = sqrt(sum(Weight.^2, 2));
amp_est = amp_est ./ sum(amp_est);
ind_est = find(amp_est > 1e-8);

amp_true = sqrt(sum(X_true.^2, 2));
amp_true = amp_true ./ sum(amp_true);

if nargin < 4 
  ind_true = find(amp_true > 1e-8);
end
  
ind_isec = union(ind_est, ind_true);
em = emd_hat_gd_metric_mex(amp_est(ind_isec), amp_true(ind_isec), D(ind_isec, ind_isec)) ./ max(D(:));

[max_co in_co] = max(abs(corr(Weight', X_true(ind_true, :)')));

mean_co = mean(max_co);

mean_co_dist = mean(diag(D(ind_true, in_co))); % in mm

%% Implement the EMD for the whole distrubution for extended sources. 



