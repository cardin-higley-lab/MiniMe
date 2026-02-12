function [dFF] = dff_bottom10(ca)

    % add lowest val to all (negative values will screw w things)
    minVal = min(ca, [], 'all');
    caAdj = ca + abs(minVal);
    
    %sort pixelwise (ie sort the frame-values per pixel)
    caSort = sort(caAdj, 2); 
    
    % get bottom 10% of vals 
    numValsInBottom10 = floor(size(caSort, 2).* 0.10); % used your bottom 10% as baseline;
    bottom10 = caSort(:, 1:numValsInBottom10);
    F0 = mean(bottom10, 2, 'omitnan') ;   % average of the bottom 10%
    dFF = (caAdj - F0)./(F0);
    dFF = single(dFF);

end  