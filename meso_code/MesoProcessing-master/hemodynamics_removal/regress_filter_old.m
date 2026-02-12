function [sig_noise1, sig_noise2, data_filt] = regress_filter(blue_sig, uv_sig, brain_mask, patchsize, inds2process)

% Regression of hempdynamics based on spatial info.
% Input: 
%  blue_sig    - calcium dependent signal - pixels X time (masked)
%  uv_sig      - calcium independent - pixels over time (masked)
%  brain_mask  - matrix of R*C = #pixels of zeors and ones indicating
%               location of brain.
% patchsize    - half size of patch + 1 for spatial info
%inds2process  - which indices to process (usefull for batch processing on farnam)
% Output: 
%   sig_noise1 - estimated standard deviation of noise in blue - pixels X 1
%   sig_noise2 - estimated standard deviation of noise in uv - pixels X 1
%   data_filt  - regressed signal
%     
% Written by:
%    Hadas Benisty, hadas.benisty@gmail.com
%    Boris Landa
%    April 6 2020
%
%------------------------------------------------------------------------













patchsize = (patchsize + 1)/2;  

% DIF CONVERT TO DOUBLE
blue_sig = double(blue_sig);
uv_sig = double(uv_sig);

[R,C] = size(brain_mask);                       % get # rows and columns of mask, typically 256x256
[Y,X] = meshgrid(1:C,1:R);                      % makes 256x256 grid, with vals = row number
X=X(:);Y=Y(:);                                  % flattens x and y, 1:256-->1:256 (256 times)
braininds = find(brain_mask);                   % gets indices of 1s in brain mask, counts down rows first

[~, blueFrames] = size(blue_sig);
[~, uvFrames] = size(uv_sig);
maxLength = min(blueFrames, uvFrames);
blue_sig = blue_sig(:, 1:maxLength);
uv_sig = uv_sig(:, 1:maxLength);

mean_blue =  mean(blue_sig,2);                      % finds the mean dff for each pixel (avgs along time dimension)
rmmean_blue = bsxfun(@minus, blue_sig, mean_blue);  % subtracts the time-averaged mean from each individual timepoint, pixelwise 

mean_uv =  mean(uv_sig,2);                          % finds the mean dff for each pixel (avgs along time dimension)
rmmean_uv = bsxfun(@minus, uv_sig, mean_uv);        % subtracts the time-averaged mean from each individual timepoint, pixelwise 

[~,Nt] = size(blue_sig);                        % Nt = columns of blue_sig, which equals number of frames (n of time)
Np = length(inds2process);                      % Np = number of pixels 
sig_noise1 = zeros(Np,1);                       % pre-allocating space for stDev for pixels (blue frames)
sig_noise2 = zeros(Np,1);                       % pre-allocating space for stDev for pixels (uv frames)
data_filt = nan(Np, Nt, 'single');

clear blue_sig
clear uv_sig

% now that everything needed is defined/preallocated, can run pixelwise correction
for piii=1:length(inds2process)    
    %piii=6479; %for testing 
    
    pii = inds2process(piii);
    if ~ismember(pii , braininds) % only continues if pixel is within brainmask
        continue;
    end

    %finds the row and column val of pii if it were in the 256x256 grid 
    x = X(pii); 
    y = Y(pii);
    
    patchPixels = find(abs(X-x) < patchsize & abs(Y-y) < patchsize);   %how it grabs the surrounding patch
    patchPixels = intersect(patchPixels, braininds); %makes sure that full patch is within the mask
    
    patch_inds=zeros(1,length(patchPixels));
    for j=1:length(patchPixels)
        patch_inds(j) = find(braininds == patchPixels(j)); %gets indices of braininds of the patch (ie where is the patch within the mask)
    end
    
%     % find if any pixels have nan at some point
%     S1 = sum(rmmean_blue(nn_inds, :),2);
%     S2 = sum(rmmean_uv(nn_inds, :),2);     
%     nn_inds=nn_inds(~isnan(S1)&~isnan(S2));    % cjb changed

    nn0_inds = find(braininds == pii); % get index of pixel in brain mask 
    
    midpixel = find(patch_inds==nn0_inds); %where is the pixel within the patch (if at the edge wont be centered)
    if isempty(midpixel)
        continue;
    end
    
    y1 = rmmean_blue(patch_inds, :); % gets pixelwise data (w time-averaged mean subtracted) for only the patch indeces 
    y2 = rmmean_uv(patch_inds, :); 

    if all(y1(:)==0)&&  all(y2(:)==0)
        continue;
    end

    d=size(y1,1); % gets the size of the patch (not always patchsize bc some pixels are at the edge)
    y = cat(1, y1, y2);
    sigy = y*y.'/size(y,2);
    sigy1 = sigy(1:d,1:d);
    sigy12 = sigy(1:d,d+1:end);
    sigy2 = sigy(1+d:end,1+d:end);
    sig_noise1(piii) = median(svd(sigy1));
    sig_noise2(piii) = median(svd(sigy2));


    sigx = sigy1-sig_noise1(piii)*eye(d) - sigy12*pinv((sigy2-sig_noise2(piii)*eye(d)))*sigy12';

    xhat = [sigx zeros(d)]*pinv(sigy)*cat(1, rmmean_blue(patch_inds, :), rmmean_uv(patch_inds, :));
    data_filt(piii,:) = xhat(midpixel,:);

end
