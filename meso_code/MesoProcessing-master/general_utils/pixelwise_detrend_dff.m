function [detrended_mat, mval, sval] = pixelwise_detrend_dff(aligned, inds2keep, filtLen)

    filtcoeff = fir1(filtLen,0.001);
    
    PixxTime_bl = aligned;
    PixxTime_bl = single(PixxTime_bl(inds2keep,:)); % subsets aligned matrix by pixels within mask
    PixxTime_bl_base = filter(filtcoeff, 1, PixxTime_bl.').'; % filters
    PixxTime_bl_base = cat(2, PixxTime_bl_base, repmat(PixxTime_bl_base(:,end), 1, filtLen/2)); % appends 75 (filtLen/2) cols to end (copies of last col)
    PixxTime_bl_base = PixxTime_bl_base(:,filtLen/2+1:end); % removes first 75 cols from baseline
    detrended = single(bsxfun(@minus, PixxTime_bl,PixxTime_bl_base)); % detrends by removing baseline
    detrended = cat(2, nan(size(detrended(:, 1:filtLen/2))), detrended(:, filtLen/2+1:end)); % converts first 75 columns into nans

    mval = mean(detrended(:, filtLen/2+1:end),2); % finds the mean of the detrended matrix
    sval = std(detrended(:, filtLen/2+1:end),[],2); % finds the standard deviation of the detrended matrix
    
    detrended_mat = single(bsxfun(@minus, detrended,mval)); % subtracts the mean from detrended mat
    detrended_mat = single(bsxfun(@rdivide, detrended_mat,sval)); % divides each trace by its standard deviation

end