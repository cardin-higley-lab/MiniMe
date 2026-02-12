function [dFF, F0, bottom10] = lowfaceNorm_bottom10(states, mesoOnTs, pixelwise_traces)

lowfaceOn = states.faceLowSitOn;
lowfaceOff = states.faceLowSitOff;
if isempty(lowfaceOn)
    error('cannot perform f0 normalization, no low facemap')
end

stateInds = false(size(pixelwise_traces,2),1);
for ii = 1:length(lowfaceOff)
    stateInds(mesoOnTs>= lowfaceOn(ii) & mesoOnTs<= lowfaceOff(ii)) = true;
end

Len = min(size(pixelwise_traces, 2), length(stateInds));
stateInds = stateInds(1:Len);
pixelwise_traces = pixelwise_traces(:, 1:Len);

lowfaceData = pixelwise_traces(:,stateInds);
lowface_sorted = sort(lowfaceData, 2);

% get bottom 10% of vals
numValsInBottom10 = floor(size(lowface_sorted, 2).* 0.10); % used your bottom 10% as baseline;
bottom10 = lowface_sorted(:, 1:numValsInBottom10);
F0 = mean(bottom10, 2, 'omitnan') ;   % average of the bottom 10%

dFF = (pixelwise_traces - F0)./F0;
dFF = single(dFF);

end