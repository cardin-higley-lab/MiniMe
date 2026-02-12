function [recon_signals, R2] = predict_behavior_from_corr(t_win, Kfolds, tC, diffmap, t_imaging, behavior_traces, J)


Len = t_imaging(t_win(end))-t_imaging(t_win(1)); 
chunksize = floor(Len/Kfolds);

% resample behavior traces
behavior_traces_resampled = cell(length(behavior_traces), 1);
for i = 1:length(behavior_traces)
    behavior_traces_resampled{i} = interp1(t_imaging(1:size(t_win,2)), behavior_traces{i}, t_imaging(t_win));
end
 

recon_signals = cell(length(behavior_traces), 1);
for n = 1:length(behavior_traces)
    recon_signals{n} = nan(size(behavior_traces_resampled{n}));
end

for chunk_i = 1:Kfolds
    %% define train and test chunks 
    teststart = t_imaging(1) + (chunk_i - 1)*chunksize; 
    testend = t_imaging(1) + chunk_i*chunksize;
    
    testinds = findClosestDouble(teststart,t_imaging(t_win)):findClosestDouble(testend,t_imaging(t_win));
    traininds = setdiff(1:length(t_win), testinds);
    
    behavior_resampled_tr = cell(length(behavior_traces_resampled),1);
    for n = 1:length(behavior_traces)
        behavior_resampled_tr{n} = behavior_traces_resampled{n}(traininds);
    end

    diffmap_tr = diffmap(:, traininds);
    diffmap_te = diffmap(:, testinds);

    %% train
    mdlglobal_phi = cell(length(behavior_traces), 1);

    for n=1:length(behavior_traces)
        mdlglobal_phi{n} = fitlm(diffmap_tr(1:J,:)', behavior_resampled_tr{n});
    end

    %% test
    isnan_te = squeeze(isnan(tC(1,1,testinds)));
   
    for n = 1:length(behavior_traces)
        recon_signals{n}(testinds(~isnan_te)) = predict(mdlglobal_phi{n}, diffmap_te(1:J,:)'); 
    end
    
    
end


R2 = max(getstats(recon_signals, behavior_traces_resampled), 0);

end


