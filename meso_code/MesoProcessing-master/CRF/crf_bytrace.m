

data_dir = 'K:\rachel\CRF\processed_08052024\CRF\';
fs = 10; preSeconds = 2; postSeconds = 5;
t = -preSeconds:(preSeconds+postSeconds)/(preSeconds*fs+postSeconds*fs):postSeconds;

response_interval = (preSeconds*fs)+1:(1*fs)+(preSeconds*fs); %for 1s post-stim onset
%response_interval = (preSeconds*fs)+1:(2*fs)+(preSeconds*fs); %for 2s post-stim onset

fig_directory = 'D:\meso_demo\figures\cdkl5_fall2024\crf\';

contrast_vals = [1 2 5 10 20 50 100];

%% load data 

% CDKL5 (+/Y) and CDKL5 (-/Y) 
control_m = { [data_dir 'ro139_365'], [data_dir 'ro159-406'], [data_dir 'ro193-604'], [data_dir 'ro230-539'], [data_dir 'ro321-277'], [data_dir 'ro337-315']}; 
mutant_m = { [data_dir 'ro139_363'],[data_dir 'ro139_366'],  [data_dir 'ro159-405'], [data_dir 'ro230-540'], [data_dir 'ro321-278'], [data_dir 'ro322-227'], [data_dir 'ro337-316']}; %

% CDKL5 (+/+) and CDKL5 (-/+)
control_f = {[data_dir 'ro159-401'], [data_dir 'ro159-403'], [data_dir 'ro159-404'], [data_dir 'ro229-535'], [data_dir 'ro229-536'], [data_dir 'ro323-231']};
het_f = {[data_dir 'ro159-402'], [data_dir 'ro229-537'], [data_dir 'ro323-229'], [data_dir 'ro323-230'], [data_dir 'ro336-311'], [data_dir 'ro336-312']}; %


% % Controls (regardless of sex) and Mutants (regardless of sex) 
% control_mf = {[data_dir 'ro139_365'], [data_dir 'ro159-406'], [data_dir 'ro193-604'], [data_dir 'ro230-539'], [data_dir 'ro159-401'], ...
%     [data_dir 'ro159-403'], [data_dir 'ro159-404'], [data_dir 'ro229-535'], [data_dir 'ro229-536']}; 
% mutant_mf = {[data_dir 'ro139_363'], [data_dir 'ro159-405'], [data_dir 'ro230-540'], [data_dir 'ro159-402'], [data_dir 'ro229-537']};


crf_struct = struct();
pre_struct = struct();
groups = {control_m, mutant_m, control_f, het_f}; %, control_mf, mutant_mf
group_names = {'control_m', 'mutant_m', 'control_f', 'het_f'}; %, 'control_mf', 'mutant_mf'
for geno = 1:length(groups)
    
    group = groups{geno}; % group = genotype/sex
    group_v1_aligned_avg_loco = nan(length(t)-1, 7, length(group));
    group_v1_aligned_avg_sit = nan(length(t)-1, 7, length(group));

    group_mean_locoResponses = nan(7, length(group));
    group_auc_locoResponses = nan(7, length(group));
    group_mean_sitResponses = nan(7, length(group));
    group_auc_sitResponses = nan(7, length(group));
    group_peak_locoResponses = nan(7, length(group));
    group_peak_sitResponses = nan(7, length(group));


    pre_struct.(group_names{geno}) = struct();
    blah=1;

    for animals = 1:length(group)
        
        animal = group(animals);
        disp(['processing: ' cell2mat(animal)])
        subfolders = dir(fullfile(cell2mat(animal)));

        animal_v1_aligned_traces_loco = [];
        animal_v1_aligned_traces_sit = [];

        for sessions = 1:length(subfolders)
            session = subfolders(sessions).name;
            if length(session) < 3
                continue;
            end
            folder = cell2mat(fullfile(animal, session));
            load(fullfile(folder, 'stim_order.mat'), 'contrasts')
            stims_per_contrast = length(contrasts)/length(unique(contrasts));

            load(fullfile(folder, 'final_timestamps.mat'), 'spike2_data');
            if isfield(spike2_data, 'blueOnTimestamps'); mesoOn = spike2_data.blueOnTimestamps; 
            else; mesoOn = spike2_data.greenOnTimestamps; end
            visOn = spike2_data.diodeOnTimestamps;
            load(fullfile(folder, 'Ca_traces_spt_patch9_lowface10p_Allen.mat'), 'parcels_time_trace');
            load(fullfile(folder, 'state_timestamps.mat'), 'states');
            
            [alignedData, ~] = averageMovieToEvent(parcels_time_trace, mesoOn, visOn, preSeconds*fs, postSeconds*fs);
            
            % zscore each trace based on m and sd of pre stimulus
            pre_traces = alignedData(:, :, 1:(preSeconds*fs)); 
            pre_traces_mean = mean(pre_traces, 3, 'omitnan');
            pre_traces_sd = std(pre_traces, 0, 3, 'omitnan');
            alignedData_z = (alignedData - pre_traces_mean)./pre_traces_sd;
           
            aligned_traces = nan(stims_per_contrast, length(t)-1, length(contrast_vals));
            for i = 1:length(contrast_vals)
                contrast_inds = contrasts==contrast_vals(i);
                if sum(contrast_inds)>0
                    aligned_traces(:,:,i) = squeeze(alignedData_z(contrast_inds, 2, :)); %change between 1 and 2 for left or right v1 
                else
                    aligned_traces(:,:,i) = nan(stims_per_contrast, length(t)-1);
                end
            end


            [~, loco] = calc_crf_responses(contrasts, parcels_time_trace, spike2_data, states);
            new_loco = reshape_loco(loco, t);
            
            v1_aligned_traces_loco = nan(size(aligned_traces));
            v1_aligned_traces_loco(new_loco==1) = aligned_traces(new_loco==1);
            animal_v1_aligned_traces_loco = cat(1, animal_v1_aligned_traces_loco, v1_aligned_traces_loco);
            
            v1_aligned_traces_sit = nan(size(aligned_traces));
            v1_aligned_traces_sit(new_loco==0) = aligned_traces(new_loco==0);
            animal_v1_aligned_traces_sit = cat(1, animal_v1_aligned_traces_sit, v1_aligned_traces_sit);

            pre_struct.(group_names{geno}).pretrace{blah} = squeeze(pre_traces(:, 1, :));
            blah = blah+1;
        end

        animal_v1_aligned_avg_loco = squeeze(mean(animal_v1_aligned_traces_loco, 1, 'omitnan'));
        animal_v1_aligned_avg_sit = squeeze(mean(animal_v1_aligned_traces_sit, 1, 'omitnan'));

        group_v1_aligned_avg_loco(:,:,animals) = animal_v1_aligned_avg_loco;
        group_v1_aligned_avg_sit(:,:,animals) = animal_v1_aligned_avg_sit;

        mean_locoResponses = nan(length(contrast_vals), 1); %auc_locoResponses = nan(length(contrast_vals), 1);
        mean_sitResponses = nan(length(contrast_vals), 1); %auc_sitResponses = nan(length(contrast_vals), 1);
        peak_locoResponses = nan(length(contrast_vals), 1); 
        peak_sitResponses = nan(length(contrast_vals), 1); 

        for i = 1:length(contrast_vals)

            trace_loco = animal_v1_aligned_avg_loco(:,i);
            mean_locoResponses(i) = mean(trace_loco(response_interval), 'omitnan');
            peak_locoResponses(i) = max(trace_loco(response_interval));
            %auc_locoResponses(i) = trapz(trace_loco(response_interval));

            trace_sit = animal_v1_aligned_avg_sit(:,i);
            mean_sitResponses(i) = mean(trace_sit(response_interval), 'omitnan');
            peak_sitResponses(i) = max(trace_sit(response_interval));
            %auc_sitResponses(i) = trapz(trace_sit(response_interval));

        end

        group_mean_locoResponses(:, animals) = mean_locoResponses; 
        group_peak_locoResponses(:, animals) = peak_locoResponses; 
        %group_auc_locoResponses(:, animals) = auc_locoResponses;
        group_mean_sitResponses(:, animals) = mean_sitResponses; 
        group_peak_sitResponses(:, animals) = peak_sitResponses; 
        %group_auc_sitResponses(:, animals) = auc_sitResponses;

    end

    crf_struct.([group_names{geno} '_locoVisTrace']) = group_v1_aligned_avg_loco;
    crf_struct.([group_names{geno} '_sitVisTrace']) = group_v1_aligned_avg_sit;
    
    crf_struct.([group_names{geno} '_locoVisResponse_mean']) = group_mean_locoResponses;
    crf_struct.([group_names{geno} '_sitVisResponse_mean']) = group_mean_sitResponses;
    %crf_struct.([group_names{geno} '_locoVisResponse_auc']) = group_auc_locoResponses;
    %crf_struct.([group_names{geno} '_sitVisResponse_auc']) = group_auc_sitResponses;
    crf_struct.([group_names{geno} '_locoVisResponse_peak']) = group_peak_locoResponses;
    crf_struct.([group_names{geno} '_sitVisResponse_peak']) = group_peak_sitResponses;

end


 











%% 100p contrast vis response - Locomotion (7=7th contrast aka 100)

control_m_mean = mean(crf_struct.control_m_locoVisTrace(:,7,:), 3, 'omitnan');
control_m_se = std(crf_struct.control_m_locoVisTrace(:,7,:), 0, 3, 'omitnan')/sqrt(size(crf_struct.control_m_locoVisTrace(:,7,:), 3));

mutant_m_mean = mean(crf_struct.mutant_m_locoVisTrace(:,7,:), 3, 'omitnan');
mutant_m_se = std(crf_struct.mutant_m_locoVisTrace(:,7,:), 0, 3, 'omitnan')/sqrt(size(crf_struct.mutant_m_locoVisTrace(:,7,:), 3));

test = figure(); hold on;
shadedErrorBar(t(2:end-20), control_m_mean(1:end-20), control_m_se(1:end-20), 'lineprops', {'k', 'markerfacecolor','k'}); 
shadedErrorBar(t(2:end-20), mutant_m_mean(1:end-20), mutant_m_se(1:end-20), 'lineprops', {'r', 'markerfacecolor','r'}); 
ylim([-2 5]); xline(0)
xlabel('Time (s)'); ylabel('v1 response'); 
title('100% contrast response  - Locomotion'); 
exportgraphics(test, [fig_directory, '\v1_100p_trace_loco_males_global_z.pdf'],'ContentType','vector')


%% 100p contrast vis response - Sitting

control_m_mean = mean(crf_struct.control_m_sitVisTrace(:,7,:), 3, 'omitnan');
control_m_se = std(crf_struct.control_m_sitVisTrace(:,7,:), 0, 3, 'omitnan')/sqrt(size(crf_struct.control_m_sitVisTrace(:,7,:), 3));

mutant_m_mean = mean(crf_struct.mutant_m_sitVisTrace(:,7,:), 3, 'omitnan');
mutant_m_se = std(crf_struct.mutant_m_sitVisTrace(:,7,:), 0, 3, 'omitnan')/sqrt(size(crf_struct.mutant_m_sitVisTrace(:,7,:), 3));

test = figure(); hold on;
shadedErrorBar(t(2:end-20), control_m_mean(1:end-20), control_m_se(1:end-20), 'lineprops', {'k', 'markerfacecolor','k'}); 
shadedErrorBar(t(2:end-20), mutant_m_mean(1:end-20), mutant_m_se(1:end-20), 'lineprops', {'r', 'markerfacecolor','r'}); 
ylim([-2 5]); xline(0)
xlabel('Time (s)'); ylabel('v1 response'); 
title('100% contrast response  - Sitting');  
exportgraphics(test, [fig_directory, '\v1_100p_trace_sitting_males_global_z.pdf'],'ContentType','vector')


%% females..  % 100p contrast vis response - Locomotion (7=7th contrast aka 100)

control_m_mean = mean(crf_struct.control_f_locoVisTrace(:,7,:), 3, 'omitnan');
control_m_se = std(crf_struct.control_f_locoVisTrace(:,7,:), 0, 3, 'omitnan')/sqrt(size(crf_struct.control_f_locoVisTrace(:,7,:), 3));

mutant_m_mean = mean(crf_struct.het_f_locoVisTrace(:,7,:), 3, 'omitnan');
mutant_m_se = std(crf_struct.het_f_locoVisTrace(:,7,:), 0, 3, 'omitnan')/sqrt(size(crf_struct.het_f_locoVisTrace(:,7,:), 3));

test = figure(); hold on;
shadedErrorBar(t(2:end-20), control_m_mean(1:end-20), control_m_se(1:end-20), 'lineprops', {'k', 'markerfacecolor','k'}); 
shadedErrorBar(t(2:end-20), mutant_m_mean(1:end-20), mutant_m_se(1:end-20), 'lineprops', {'r', 'markerfacecolor','r'}); 
ylim([-2 5]); xline(0)
xlabel('Time (s)'); ylabel('v1 response'); 
title('100% contrast response  - Locomotion - Females'); 
exportgraphics(test, [fig_directory, '\v1_100p_trace_loco_females_global_z.pdf'],'ContentType','vector')


%%  females..  100p contrast vis response - Sitting

control_m_mean = mean(crf_struct.control_f_sitVisTrace(:,7,:), 3, 'omitnan');
control_m_se = std(crf_struct.control_f_sitVisTrace(:,7,:), 0, 3, 'omitnan')/sqrt(size(crf_struct.control_m_sitVisTrace(:,7,:), 3));

mutant_m_mean = mean(crf_struct.het_f_sitVisTrace(:,7,:), 3, 'omitnan');
mutant_m_se = std(crf_struct.het_f_sitVisTrace(:,7,:), 0, 3, 'omitnan')/sqrt(size(crf_struct.mutant_m_sitVisTrace(:,7,:), 3));

test = figure(); hold on;
shadedErrorBar(t(2:end-20), control_m_mean(1:end-20), control_m_se(1:end-20), 'lineprops', {'k', 'markerfacecolor','k'}); 
shadedErrorBar(t(2:end-20), mutant_m_mean(1:end-20), mutant_m_se(1:end-20), 'lineprops', {'r', 'markerfacecolor','r'}); 
ylim([-2 5]); xline(0)
xlabel('Time (s)'); ylabel('v1 response'); 
title('100% contrast response  - Sitting - Females');  
exportgraphics(test, [fig_directory, '\v1_100p_trace_sit_females_global_z.pdf'],'ContentType','vector')














%% CRFs - males - loco

control_m_mean = mean(crf_struct.control_m_locoVisResponse_mean, 2, 'omitnan');
control_m_se = std(crf_struct.control_m_locoVisResponse_mean, 0, 2, 'omitnan')/sqrt(size(crf_struct.control_m_locoVisResponse_mean, 2));

mutant_m_mean = mean(crf_struct.mutant_m_locoVisResponse_mean, 2, 'omitnan');
mutant_m_se = std(crf_struct.mutant_m_locoVisResponse_mean, 0, 2, 'omitnan')/sqrt(size(crf_struct.mutant_m_locoVisResponse_mean, 2));

contrasts = ([1, 2, 5, 10, 20, 50, 100]);
test = figure(); hold on;
errorbar(contrasts, control_m_mean, control_m_se, 'k'); plot(contrasts, control_m_mean, 'k.', 'MarkerSize', 12)
errorbar(contrasts, mutant_m_mean, mutant_m_se, 'r'); plot(contrasts, mutant_m_mean, 'r.', 'MarkerSize', 12)
ylim([-1 2]); xlabel('Contrast ( %)'); ylabel('mean v1 response'); 
title('CRF - Locomotion - Males'); 
exportgraphics(test, [fig_directory, '\v1_crf_loco_males.pdf'],'ContentType','vector')




%% CRFs - males - sit

control_m_mean = mean(crf_struct.control_m_sitVisResponse_mean, 2, 'omitnan');
control_m_se = std(crf_struct.control_m_sitVisResponse_mean, 0, 2, 'omitnan')/sqrt(size(crf_struct.control_m_sitVisResponse_mean, 2));

mutant_m_mean = mean(crf_struct.mutant_m_sitVisResponse_mean, 2, 'omitnan');
mutant_m_se = std(crf_struct.mutant_m_sitVisResponse_mean, 0, 2, 'omitnan')/sqrt(size(crf_struct.mutant_m_sitVisResponse_mean, 2));

contrasts = ([1, 2, 5, 10, 20, 50, 100]);
test = figure(); hold on;
errorbar(contrasts, control_m_mean, control_m_se, 'k'); plot(contrasts, control_m_mean, 'k.', 'MarkerSize', 12)
errorbar(contrasts, mutant_m_mean, mutant_m_se, 'r'); plot(contrasts, mutant_m_mean, 'r.', 'MarkerSize', 12)
ylim([-1 2]); xlabel('Contrast (%)'); ylabel('mean V1 response'); 
title('CRF - Not Loco - Males'); 
exportgraphics(test, [fig_directory, 'v1_crf_sit_males.pdf'],'ContentType','vector')






%% CRFs - females - loco

control_f_mean = mean(crf_struct.control_f_locoVisResponse_mean, 2, 'omitnan');
control_f_se = std(crf_struct.control_f_locoVisResponse_mean, 0, 2, 'omitnan')/sqrt(size(crf_struct.control_f_locoVisResponse_mean, 2));

mutant_f_mean = mean(crf_struct.het_f_locoVisResponse_mean, 2, 'omitnan');
mutant_f_se = std(crf_struct.het_f_locoVisResponse_mean, 0, 2, 'omitnan')/sqrt(size(crf_struct.het_f_locoVisResponse_mean, 2));

contrasts = ([1, 2, 5, 10, 20, 50, 100]);
test = figure(); hold on;
errorbar(contrasts, control_f_mean, control_f_se, 'k'); plot(contrasts, control_f_mean, 'k.', 'MarkerSize', 12)
errorbar(contrasts, mutant_f_mean, mutant_f_se, 'r'); plot(contrasts, mutant_f_mean, 'r.', 'MarkerSize', 12)
ylim([-1 2]); xlabel('Contrast ( %)'); ylabel('mean v1 response'); 
title('CRF - Locomotion - Females'); 
exportgraphics(test, [fig_directory, '\v1_crf_loco_females.pdf'],'ContentType','vector')


%% CRFs - females - sit

control_f_mean = mean(crf_struct.control_f_sitVisResponse_mean, 2, 'omitnan');
control_f_se = std(crf_struct.control_f_sitVisResponse_mean, 0, 2, 'omitnan')/sqrt(size(crf_struct.control_f_sitVisResponse_mean, 2));

mutant_f_mean = mean(crf_struct.het_f_sitVisResponse_mean, 2, 'omitnan');
mutant_f_se = std(crf_struct.het_f_sitVisResponse_mean, 0, 2, 'omitnan')/sqrt(size(crf_struct.het_f_sitVisResponse_mean, 2));

contrasts = ([1, 2, 5, 10, 20, 50, 100]);
test = figure(); hold on;
errorbar(contrasts, control_f_mean, control_f_se, 'k'); plot(contrasts, control_f_mean, 'k.', 'MarkerSize', 12)
errorbar(contrasts, mutant_f_mean, mutant_f_se, 'r'); plot(contrasts, mutant_f_mean, 'r.', 'MarkerSize', 12)
ylim([-1 2]); xlabel('Contrast ( %)'); ylabel('mean v1 response'); 
title('CRF - Not Sit - Females'); 
exportgraphics(test, [fig_directory, '\v1_crf_sit_females.pdf'],'ContentType','vector')




%%


control_m = cat(1, pre_struct.control_m.pretrace{:});
mutant_m = cat(1, pre_struct.mutant_m.pretrace{:});

control_f = cat(1, pre_struct.control_f.pretrace{:});
het_f = cat(1, pre_struct.het_f.pretrace{:});

figure();histogram(control_m(:)); title('control_m'); xlim([-0.05 0.5])
figure();histogram(mutant_m(:)); title('mutant_m'); xlim([-0.05 0.5])
figure();histogram(control_f(:)); title('control_f'); xlim([-0.05 0.5])
figure();histogram(het_f(:)); title('het_f'); xlim([-0.05 0.5])

