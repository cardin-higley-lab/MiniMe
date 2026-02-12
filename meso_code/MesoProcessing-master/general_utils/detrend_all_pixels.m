function [outputsig, baseline] = detrend_all_pixels(insig, paramsDeterend)

filtLen = paramsDeterend.filtLen;
filtcutoff = paramsDeterend.filtcutoff;
method = paramsDeterend.method;


if nargin == 1
    filtLen =1000;filtcutoff = 0.001;
    method = 'FIR';
end
switch method
    case 'TOPHAT'
        error('Needs some work');
    case 'FIR'        
        filtcoeff = fir1(filtLen,filtcutoff);
        baseline = filter(filtcoeff, 1, insig.').';    
        baseline = baseline(:, filtLen:end);%remove 1/2 filtlen to account for temporal delay and another 1/2 filtlen to account for edge artifact 
        insig = insig(:, filtLen/2:(end-(filtLen/2)));%remove the first 1/2 filteln frames to account for edge artifact and the last 1/2 filtlen artifacts because of temporal delay in baseline 
        outputsig=insig - baseline;
        %outputsig=bsxfun(@minus, outputsig, min(outputsig));%% make values positive 
    case 'SimpleMean'
        baseline = ones(size(insig)).*mean(insig,1);
        outputsig = insig-baseline;
end



