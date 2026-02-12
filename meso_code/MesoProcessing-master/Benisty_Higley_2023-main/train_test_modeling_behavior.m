%% main function for modeling behavior 
% inputs: 
%   - imaging_signals - a matrix of parcels over time
%   - t_imaging       - time stamps for the imaging signal
%   - behavior_traces - a matrix of behavior signals over time
%   - phi_c           - a matrix of diff map elements over time windows (diffmap)
%   - t_win           - a vector of mid points of time windows corresponding to phi_c


%  res_xraw_phi_ridge, res_xrawsh_phi_ridge and res_xraw_phish_ridge each
%  have 4 fields: 
%      R_tr = R2 for train
%      R_te = R2 for test
%      C_tr = correlation between real behavior and estimated for train
%      C_te = correlation between real behavior and estimated for test

function [res_xraw_phi, res_xrawsh_phi, res_xraw_phish, predicted_xraw_phi_ridge, predicted_xrawsh_phi_ridge, predicted_xraw_phish_ridge] = train_test_modeling_behavior(imaging_signals, t_imaging, behavior_traces, phi_c, t_win, J, Kfolds)

t_imaging = t_imaging(1:size(imaging_signals,2));  
phi_c = phi_c(:, t_win<=length(t_imaging));
t_win = t_win(t_win<=length(t_imaging));

% resample behavior traces
behavior_traces_resampled = cell(length(behavior_traces), 1);
for i = 1:length(behavior_traces)
    behavior_traces_resampled{i} = interp1(t_imaging, behavior_traces{i}, t_imaging(t_win));
end
behavior_traces_resampled = cell2mat(behavior_traces_resampled');

% resample imaging signals activity to match mid points of windows of phi_c
imaging_signalsc = interp1(t_imaging, imaging_signals', t_imaging(t_win))';

% concatenating activity and phi_c
X = cat(1, zscore(imaging_signalsc')', phi_c(1:J,:)); 

% indices of activity in X
xinds = 1:size(imaging_signalsc, 1);

% ridge regression of X + phi_c
[res_xraw_phi, ~, ~, predicted_xraw_phi_ridge] = train_test_ridge(X, behavior_traces_resampled, Kfolds);

% ridge regression of X(shuffled) + phi_c
[res_xrawsh_phi, ~, ~, predicted_xrawsh_phi_ridge] = train_test_ridge(X, behavior_traces_resampled, Kfolds, xinds);

% ridge regression of X + phi_c(shuffled)
[res_xraw_phish, ~, ~, predicted_xraw_phish_ridge] = train_test_ridge(X, behavior_traces_resampled, Kfolds, setdiff(1:size(X, 1), xinds));

end




