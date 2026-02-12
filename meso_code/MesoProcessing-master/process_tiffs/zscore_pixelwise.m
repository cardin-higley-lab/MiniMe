function [z_blue, mval, sval] = zscore_pixelwise(blue, filtLen)


mval = mean(blue(:, filtLen/2+1:end),2, 'omitnan');
sval = std(blue(:, filtLen/2+1:end),[],2, 'omitnan');

z_blue = single(bsxfun(@minus, blue, mval));
z_blue = single(bsxfun(@rdivide, z_blue, sval));
