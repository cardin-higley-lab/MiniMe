function [tC, t_win] = sliding_window_corr_reg(imaging_data, winsizeSec, winhopSec)

    % add tau as the noise level
    T = winsizeSec*imaging_data.fsample;
    s = svd(imaging_data.data(:,1:round(T))); 
    th = quantile(s,0.3);
    tau1 = median(s(s<=th));
    [tC, t_win] = sliding_window_corr(imaging_data, winsizeSec, winhopSec, 2*tau1);

end

