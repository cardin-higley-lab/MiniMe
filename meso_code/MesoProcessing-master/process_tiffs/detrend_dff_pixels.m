function det_sig = detrend_dff_pixels(sig_raw, paramsdeterend)

filtcoeff = fir1(paramsdeterend.filtLen,paramsdeterend.filtcutoff);

PixxTime_bl_base = filter(filtcoeff, 1, sig_raw.').';
PixxTime_bl_base = cat(2, PixxTime_bl_base, repmat(PixxTime_bl_base(:,end), 1, paramsdeterend.filtLen/2));
PixxTime_bl_base = PixxTime_bl_base(:,paramsdeterend.filtLen/2+1:end);
det_sig = single(bsxfun(@minus, sig_raw,PixxTime_bl_base)./PixxTime_bl_base);

