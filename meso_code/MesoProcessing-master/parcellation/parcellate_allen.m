function parcels_time_trace = parcellate_allen(parcells_new, normed_sig)

parcels_time_trace = zeros(length(parcells_new.names), size(normed_sig,2));
for par_i = 1:length(parcells_new.names)

    roiinds = parcells_new.indicators(:,:,par_i)==1;
    pardata = normed_sig(roiinds, :);
    detect = isoutlier(pardata);
    keep_data = nan(size(pardata));
    keep_data(detect==0) = pardata(detect==0);
    parcels_time_trace(par_i,:) = mean(keep_data, 'omitnan');

end


end