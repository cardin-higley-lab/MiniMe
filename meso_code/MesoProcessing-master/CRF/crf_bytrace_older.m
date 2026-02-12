

data_dir = 'K:\rachel\CRF\processed_08052024\CRF\';
fs = 10; preSeconds = 1; postSeconds = 5;
t = -preSeconds:(preSeconds+postSeconds)/(preSeconds*fs+postSeconds*fs):postSeconds;

pre_interval = 1:(preSeconds*fs);
response_interval = (preSeconds*fs)+1:(1*fs)+(preSeconds*fs); %for 1s post-stim onset
%response_interval = (preSeconds*fs)+1:(2*fs)+(preSeconds*fs); %for 2s post-stim onset

fig_directory = 'D:\meso_demo\figures\cdkl5_poster\';

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
            load(fullfile(folder, 'final_timestamps.mat'), 'spike2_data');
            if isfield(spike2_data, 'blueOnTimestamps'); mesoOn = spike2_data.blueOnTimestamps; 
            else; mesoOn = spike2_data.greenOnTimestamps; end
            visOn = spike2_data.diodeOnTimestamps;
            load(fullfile(folder, 'Ca_traces_spt_patch9_lowface10p_Allen.mat'), 'parcels_time_trace');
            load(fullfile(folder, 'state_timestamps.mat'), 'states');

            [~, loco] = calc_crf_responses(contrasts, parcels_time_trace, spike2_data, states);
            new_loco = reshape_loco(loco, t);
            [alignedData, preData] = averageMovieToEvent(parcels_time_trace, mesoOn, visOn, preSeconds*fs, postSeconds*fs);
            normed_aligned = alignedData - mean(preData, 3); %decide if subtract or divide... 
            detect_outliers = isoutlier(normed_aligned, "mean"); 
            sum_nan_vals = squeeze(sum(detect_outliers, 2:3));
            outlier_stims = find(sum_nan_vals > 300);
            normed_aligned(outlier_stims, :, :) = nan;

            stims_per_contrast = length(contrasts)/length(unique(contrasts));
            aligned_traces = nan(stims_per_contrast, length(t)-1, length(contrast_vals));
            for i = 1:length(contrast_vals)
                contrast_inds = contrasts==contrast_vals(i);
                if sum(contrast_inds)>0
                    aligned_traces(:,:,i) = squeeze(normed_aligned(contrast_inds, 2, :)); %change between 1 and 2 for left or right v1 
                else
                    aligned_traces(:,:,i) = nan(stims_per_contrast, length(t)-1);
                end
            end

            v1_aligned_traces_loco = nan(size(aligned_traces));
            v1_aligned_traces_loco(new_loco==1) = aligned_traces(new_loco==1);
            animal_v1_aligned_traces_loco = cat(1, animal_v1_aligned_traces_loco, v1_aligned_traces_loco);
            
            v1_aligned_traces_sit = nan(size(aligned_traces));
            v1_aligned_traces_sit(new_loco==0) = aligned_traces(new_loco==0);
            animal_v1_aligned_traces_sit = cat(1, animal_v1_aligned_traces_sit, v1_aligned_traces_sit);

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

%% visualize some of the traces...

% figure(); hold on
% plot(t(2:end), crf_struct.control_m_locoVisTrace(:,7,1), 'k'); xline(0)
% plot(t(2:end), crf_struct.control_m_locoVisTrace(:,7,2), 'k'); 
% plot(t(2:end), crf_struct.control_m_locoVisTrace(:,7,3), 'k'); 
% plot(t(2:end), crf_struct.control_m_locoVisTrace(:,7,4), 'k'); 
% 
% plot(t(2:end), crf_struct.mutant_m_locoVisTrace(:,7,1), 'r'); 
% plot(t(2:end), crf_struct.mutant_m_locoVisTrace(:,7,2), 'r'); 
% plot(t(2:end), crf_struct.mutant_m_locoVisTrace(:,7,3), 'r'); 
% 
% figure(); hold on
% plot(t(2:end), crf_struct.control_m_sitVisTrace(:,7,1), 'k'); xline(0)
% plot(t(2:end), crf_struct.control_m_sitVisTrace(:,7,2), 'k'); 
% plot(t(2:end), crf_struct.control_m_sitVisTrace(:,7,3), 'k'); 
% plot(t(2:end), crf_struct.control_m_sitVisTrace(:,7,4), 'k'); 
% 
% plot(t(2:end), crf_struct.mutant_m_sitVisTrace(:,7,1), 'r'); 
% plot(t(2:end), crf_struct.mutant_m_sitVisTrace(:,7,2), 'r'); 
% plot(t(2:end), crf_struct.mutant_m_sitVisTrace(:,7,3), 'r'); 

 
%% calculate vis response (mean) - Locomotion
 
control_m_mean = mean(crf_struct.control_m_locoVisResponse_mean, 2, 'omitnan');
control_m_se = std(crf_struct.control_m_locoVisResponse_mean, 0, 2, 'omitnan')/sqrt(size(crf_struct.control_m_locoVisResponse_mean, 2));

mutant_m_mean = mean(crf_struct.mutant_m_locoVisResponse_mean, 2, 'omitnan');
mutant_m_se = std(crf_struct.mutant_m_locoVisResponse_mean, 0, 2, 'omitnan')/sqrt(size(crf_struct.mutant_m_locoVisResponse_mean, 2));

contrasts = ([1, 2, 5, 10, 20, 50, 100]);
test = figure(); hold on;
shadedErrorBar(contrasts, control_m_mean, control_m_se, 'lineprops', {'k', 'markerfacecolor','k'}); 
plot(contrasts, control_m_mean, 'ko')
shadedErrorBar(contrasts, mutant_m_mean, mutant_m_se, 'lineprops', {'r', 'markerfacecolor','r'}); 
plot(contrasts, mutant_m_mean, 'ro');  ylim([-0.003 0.012]) %ylim([0.9 1.4]) %
xlabel('Contrast (%)'); ylabel('mean V1 response'); 
title('CRF - Locomotion - Males'); 
%exportgraphics(test, [fig_directory, '\crf_cdkl5_m_loco_1s_v1.pdf'],'ContentType','vector')


%% calculate vis response (mean) - Sitting

control_m_mean = mean(crf_struct.control_m_sitVisResponse_mean, 2, 'omitnan');
control_m_se = std(crf_struct.control_m_sitVisResponse_mean, 0, 2, 'omitnan')/sqrt(size(crf_struct.control_m_sitVisResponse_mean, 2));

mutant_m_mean = mean(crf_struct.mutant_m_sitVisResponse_mean, 2, 'omitnan');
mutant_m_se = std(crf_struct.mutant_m_sitVisResponse_mean, 0, 2, 'omitnan')/sqrt(size(crf_struct.mutant_m_sitVisResponse_mean, 2));

contrasts = ([1, 2, 5, 10, 20, 50, 100]);
test = figure(); hold on;
shadedErrorBar(contrasts, control_m_mean, control_m_se, 'lineprops', {'k', 'markerfacecolor','k'}); 
plot(contrasts, control_m_mean, 'ko')
shadedErrorBar(contrasts, mutant_m_mean, mutant_m_se, 'lineprops', {'r', 'markerfacecolor','r'}); 
plot(contrasts, mutant_m_mean, 'ro'); ylim([-0.003 0.012]) %ylim([0.9 1.4]) %ylim([-0.004 0.008])
xlabel('Contrast (%)'); ylabel('mean V1 response'); 
title('CRF - Not Loco - Males'); 
%exportgraphics(test, [fig_directory, 'crf_cdkl5_m_sitting_1s_V1.pdf'],'ContentType','vector')



%% calculate vis response (mean) - loco FEMALES - 

control_f_mean = mean(crf_struct.control_f_locoVisResponse_mean, 2, 'omitnan');
control_f_se = std(crf_struct.control_f_locoVisResponse_mean, 0, 2, 'omitnan')/sqrt(size(crf_struct.control_f_locoVisResponse_mean, 2));

mutant_f_mean = mean(crf_struct.het_f_locoVisResponse_mean, 2, 'omitnan');
mutant_f_se = std(crf_struct.het_f_locoVisResponse_mean, 0, 2, 'omitnan')/sqrt(size(crf_struct.het_f_locoVisResponse_mean, 2));

contrasts = ([1, 2, 5, 10, 20, 50, 100]);
test = figure(); hold on;
shadedErrorBar(contrasts, control_f_mean, control_f_se, 'lineprops', {'k', 'markerfacecolor','k'}); 
plot(contrasts, control_f_mean, 'ko')
shadedErrorBar(contrasts, mutant_f_mean, mutant_f_se, 'lineprops', {'r', 'markerfacecolor','r'}); 
plot(contrasts, mutant_f_mean, 'ro'); ylim([-0.003 0.012]) %ylim([0.9 1.6]) %
xlabel('Contrast (%)'); ylabel('mean V1 response'); 
title('CRF - Locomotion - Females'); 
%exportgraphics(test, 'D:\meso_demo\figures\crfoutput_09292023\crf_cdkl5_f_loco_1s_VISal.pdf','ContentType','vector')



%% calculate vis response (mean) - Sitting FEMALES - 

control_f_mean = mean(crf_struct.control_f_sitVisResponse_mean, 2, 'omitnan');
control_f_se = std(crf_struct.control_f_sitVisResponse_mean, 0, 2, 'omitnan')/sqrt(size(crf_struct.control_f_sitVisResponse_mean, 2));

mutant_f_mean = mean(crf_struct.het_f_sitVisResponse_mean, 2, 'omitnan');
mutant_f_se = std(crf_struct.het_f_sitVisResponse_mean, 0, 2, 'omitnan')/sqrt(size(crf_struct.het_f_sitVisResponse_mean, 2));

contrasts = ([1, 2, 5, 10, 20, 50, 100]);
test = figure(); hold on;
shadedErrorBar(contrasts, control_f_mean, control_f_se, 'lineprops', {'k', 'markerfacecolor','k'}); 
plot(contrasts, control_f_mean, 'ko')
shadedErrorBar(contrasts, mutant_f_mean, mutant_f_se, 'lineprops', {'r', 'markerfacecolor','r'}); 
plot(contrasts, mutant_f_mean, 'ro'); ylim([-0.003 0.012]) %ylim([0.9 1.6]) %
xlabel('Contrast (%)'); ylabel('mean V1 response'); 
title('CRF - Not Loco - Females'); 
%exportgraphics(test, 'D:\meso_demo\figures\crfoutput_09292023\crf_cdkl5_f_sitting_1s_VISal.pdf','ContentType','vector')










%% calculate vis response (peak) - Locomotion

control_m_mean = mean(crf_struct.control_m_locoVisResponse_peak, 2, 'omitnan');
control_m_se = std(crf_struct.control_m_locoVisResponse_peak, 0, 2, 'omitnan')/sqrt(size(crf_struct.control_m_locoVisResponse_peak, 2));

mutant_m_mean = mean(crf_struct.mutant_m_locoVisResponse_peak, 2, 'omitnan');
mutant_m_se = std(crf_struct.mutant_m_locoVisResponse_peak, 0, 2, 'omitnan')/sqrt(size(crf_struct.mutant_m_locoVisResponse_peak, 2));

contrasts = ([1, 2, 5, 10, 20, 50, 100]);
test = figure(); hold on;
shadedErrorBar(contrasts, control_m_mean, control_m_se, 'lineprops', {'k', 'markerfacecolor','k'}); 
plot(contrasts, control_m_mean, 'ko')
shadedErrorBar(contrasts, mutant_m_mean, mutant_m_se, 'lineprops', {'r', 'markerfacecolor','r'}); 
plot(contrasts, mutant_m_mean, 'ro');  ylim([-0.002 0.015]) %ylim([0.9 1.4]) %
xlabel('Contrast (%)'); ylabel('peak v1 response'); 
title('CRF - Locomotion - males'); 
%exportgraphics(test, 'D:\meso_demo\figures\crfoutput_09292023\crf_cdkl5_m_loco_bldiv.pdf','ContentType','vector')


%% calculate vis response (peak) - Sitting

control_m_mean = mean(crf_struct.control_m_sitVisResponse_peak, 2, 'omitnan');
control_m_se = std(crf_struct.control_m_sitVisResponse_peak, 0, 2, 'omitnan')/sqrt(size(crf_struct.control_m_sitVisResponse_peak, 2));

mutant_m_mean = mean(crf_struct.mutant_m_sitVisResponse_peak, 2, 'omitnan');
mutant_m_se = std(crf_struct.mutant_m_sitVisResponse_peak, 0, 2, 'omitnan')/sqrt(size(crf_struct.mutant_m_sitVisResponse_peak, 2));

contrasts = ([1, 2, 5, 10, 20, 50, 100]);
test = figure(); hold on;
shadedErrorBar(contrasts, control_m_mean, control_m_se, 'lineprops', {'k', 'markerfacecolor','k'}); 
plot(contrasts, control_m_mean, 'ko')
shadedErrorBar(contrasts, mutant_m_mean, mutant_m_se, 'lineprops', {'r', 'markerfacecolor','r'}); 
plot(contrasts, mutant_m_mean, 'ro'); ylim([-0.002 0.015]) %ylim([0.9 1.4]) %ylim([-0.004 0.008])
xlabel('Contrast (%)'); ylabel('peak v1 response'); 
title('CRF - Not Loco'); 
%exportgraphics(test, 'D:\meso_demo\figures\crfoutput_09292023\crf_cdkl5_m_sitting_bldiv.pdf','ContentType','vector')








%%


%% 100p contrast vis response - Locomotion (7=7th contrast aka 100)

control_m_mean = mean(crf_struct.control_m_locoVisTrace(:,7,:), 3, 'omitnan');
control_m_se = std(crf_struct.control_m_locoVisTrace(:,7,:), 0, 3, 'omitnan')/sqrt(size(crf_struct.control_m_locoVisTrace(:,7,:), 3));

mutant_m_mean = mean(crf_struct.mutant_m_locoVisTrace(:,7,:), 3, 'omitnan');
mutant_m_se = std(crf_struct.mutant_m_locoVisTrace(:,7,:), 0, 3, 'omitnan')/sqrt(size(crf_struct.mutant_m_locoVisTrace(:,7,:), 3));

test = figure(); hold on;
shadedErrorBar(t(2:end-20), control_m_mean(1:end-20), control_m_se(1:end-20), 'lineprops', {'k', 'markerfacecolor','k'}); 
shadedErrorBar(t(2:end-20), mutant_m_mean(1:end-20), mutant_m_se(1:end-20), 'lineprops', {'r', 'markerfacecolor','r'}); 
ylim([-0.01 0.015]); xline(0)
xlabel('Time (s)'); ylabel('v1 response'); 
title('100% contrast response  - Locomotion'); 
%exportgraphics(test, [fig_directory, '\cdkl5_100p_timetrace_loco_v1_male.pdf'],'ContentType','vector')


%% 100p contrast vis response - Sitting

control_m_mean = mean(crf_struct.control_m_sitVisTrace(:,7,:), 3, 'omitnan');
control_m_se = std(crf_struct.control_m_sitVisTrace(:,7,:), 0, 3, 'omitnan')/sqrt(size(crf_struct.control_m_sitVisTrace(:,7,:), 3));

mutant_m_mean = mean(crf_struct.mutant_m_sitVisTrace(:,7,:), 3, 'omitnan');
mutant_m_se = std(crf_struct.mutant_m_sitVisTrace(:,7,:), 0, 3, 'omitnan')/sqrt(size(crf_struct.mutant_m_sitVisTrace(:,7,:), 3));

test = figure(); hold on;
shadedErrorBar(t(2:end-20), control_m_mean(1:end-20), control_m_se(1:end-20), 'lineprops', {'k', 'markerfacecolor','k'}); 
shadedErrorBar(t(2:end-20), mutant_m_mean(1:end-20), mutant_m_se(1:end-20), 'lineprops', {'r', 'markerfacecolor','r'}); 
ylim([-0.01 0.015]); xline(0)
xlabel('Time (s)'); ylabel('v1 response'); 
title('100% contrast response  - Sitting');  
%exportgraphics(test, [fig_directory, '\cdkl5_100p_timetrace_sitting_v1_male.pdf'],'ContentType','vector')

%% females..  % 100p contrast vis response - Locomotion (7=7th contrast aka 100)

control_m_mean = mean(crf_struct.control_f_locoVisTrace(:,7,:), 3, 'omitnan');
control_m_se = std(crf_struct.control_f_locoVisTrace(:,7,:), 0, 3, 'omitnan')/sqrt(size(crf_struct.control_f_locoVisTrace(:,7,:), 3));

mutant_m_mean = mean(crf_struct.het_f_locoVisTrace(:,7,:), 3, 'omitnan');
mutant_m_se = std(crf_struct.het_f_locoVisTrace(:,7,:), 0, 3, 'omitnan')/sqrt(size(crf_struct.het_f_locoVisTrace(:,7,:), 3));

test = figure(); hold on;
shadedErrorBar(t(2:end-20), control_m_mean(1:end-20), control_m_se(1:end-20), 'lineprops', {'k', 'markerfacecolor','k'}); 
shadedErrorBar(t(2:end-20), mutant_m_mean(1:end-20), mutant_m_se(1:end-20), 'lineprops', {'r', 'markerfacecolor','r'}); 
ylim([-0.01 0.02]); xline(0)
xlabel('Time (s)'); ylabel('v1 response'); 
title('100% contrast response  - Locomotion - Females'); 
%exportgraphics(test, [fig_directory, '\cdkl5_100p_timetrace_loco_v1_male.pdf'],'ContentType','vector')


%%  females..  100p contrast vis response - Sitting

control_m_mean = mean(crf_struct.control_f_sitVisTrace(:,7,:), 3, 'omitnan');
control_m_se = std(crf_struct.control_f_sitVisTrace(:,7,:), 0, 3, 'omitnan')/sqrt(size(crf_struct.control_m_sitVisTrace(:,7,:), 3));

mutant_m_mean = mean(crf_struct.het_f_sitVisTrace(:,7,:), 3, 'omitnan');
mutant_m_se = std(crf_struct.het_f_sitVisTrace(:,7,:), 0, 3, 'omitnan')/sqrt(size(crf_struct.mutant_m_sitVisTrace(:,7,:), 3));

test = figure(); hold on;
shadedErrorBar(t(2:end-20), control_m_mean(1:end-20), control_m_se(1:end-20), 'lineprops', {'k', 'markerfacecolor','k'}); 
shadedErrorBar(t(2:end-20), mutant_m_mean(1:end-20), mutant_m_se(1:end-20), 'lineprops', {'r', 'markerfacecolor','r'}); 
ylim([-0.01 0.02]); xline(0)
xlabel('Time (s)'); ylabel('v1 response'); 
title('100% contrast response  - Sitting - Females');  
%exportgraphics(test, [fig_directory, '\cdkl5_100p_timetrace_sitting_v1_male.pdf'],'ContentType','vector')














%% trying something else for showing CRFs

control_m_mean = mean(crf_struct.control_m_locoVisResponse_mean, 2, 'omitnan');
control_m_se = std(crf_struct.control_m_locoVisResponse_mean, 0, 2, 'omitnan')/sqrt(size(crf_struct.control_m_locoVisResponse_mean, 2));

mutant_m_mean = mean(crf_struct.mutant_m_locoVisResponse_mean, 2, 'omitnan');
mutant_m_se = std(crf_struct.mutant_m_locoVisResponse_mean, 0, 2, 'omitnan')/sqrt(size(crf_struct.mutant_m_locoVisResponse_mean, 2));

contrasts = ([1, 2, 5, 10, 20, 50, 100]);
test = figure(); hold on;
errorbar(contrasts, control_m_mean, control_m_se, 'k'); plot(contrasts, control_m_mean, 'k.', 'MarkerSize', 12)
errorbar(contrasts, mutant_m_mean, mutant_m_se, 'r'); plot(contrasts, mutant_m_mean, 'r.', 'MarkerSize', 12)
ylim([-0.003 0.012]) ; xlabel('Contrast ( %)'); ylabel('mean v1 response'); 
title('CRF - Locomotion - Males'); 
exportgraphics(test, [fig_directory, '\crf_cdkl5_m_loco_1s_v1_newnew.pdf'],'ContentType','vector')




%% calculate vis response (mean) - Sitting

control_m_mean = mean(crf_struct.control_m_sitVisResponse_mean, 2, 'omitnan');
control_m_se = std(crf_struct.control_m_sitVisResponse_mean, 0, 2, 'omitnan')/sqrt(size(crf_struct.control_m_sitVisResponse_mean, 2));

mutant_m_mean = mean(crf_struct.mutant_m_sitVisResponse_mean, 2, 'omitnan');
mutant_m_se = std(crf_struct.mutant_m_sitVisResponse_mean, 0, 2, 'omitnan')/sqrt(size(crf_struct.mutant_m_sitVisResponse_mean, 2));

contrasts = ([1, 2, 5, 10, 20, 50, 100]);
test = figure(); hold on;
errorbar(contrasts, control_m_mean, control_m_se, 'k'); plot(contrasts, control_m_mean, 'k.', 'MarkerSize', 12)
errorbar(contrasts, mutant_m_mean, mutant_m_se, 'r'); plot(contrasts, mutant_m_mean, 'r.', 'MarkerSize', 12)
ylim([-0.003 0.012]); xlabel('Contrast (%)'); ylabel('mean V1 response'); 
title('CRF - Not Loco - Males'); 
exportgraphics(test, [fig_directory, 'crf_cdkl5_m_sitting_1s_V1_newnew.pdf'],'ContentType','vector')






%% trying something else for showing CRFs

control_m_mean = mean(crf_struct.control_m_locoVisResponse_mean, 2, 'omitnan');
control_m_se = std(crf_struct.control_m_locoVisResponse_mean, 0, 2, 'omitnan')/sqrt(size(crf_struct.control_m_locoVisResponse_mean, 2));

mutant_m_mean = mean(crf_struct.mutant_m_locoVisResponse_mean, 2, 'omitnan');
mutant_m_se = std(crf_struct.mutant_m_locoVisResponse_mean, 0, 2, 'omitnan')/sqrt(size(crf_struct.mutant_m_locoVisResponse_mean, 2));

contrasts = ([1, 2, 5, 10, 20, 50, 100]);
test = figure(); hold on;
errorbar(contrasts, control_m_mean, control_m_se, 'k'); plot(contrasts, control_m_mean, 'k.', 'MarkerSize', 12)
errorbar(contrasts, mutant_m_mean, mutant_m_se, 'r'); plot(contrasts, mutant_m_mean, 'r.', 'MarkerSize', 12)
ylim([-0.003 0.012]) ; xlabel('Contrast ( %)'); ylabel('mean v1 response'); 
title('CRF - Locomotion - Males'); 
exportgraphics(test, [fig_directory, '\crf_cdkl5_m_loco_1s_v1_newnew.pdf'],'ContentType','vector')



%% combined males and female plot


%% calculate vis response (mean) - Sitting

control_m_mean = mean(crf_struct.control_m_sitVisResponse_mean, 2, 'omitnan');
control_m_se = std(crf_struct.control_m_sitVisResponse_mean, 0, 2, 'omitnan')/sqrt(size(crf_struct.control_m_sitVisResponse_mean, 2));

mutant_m_mean = mean(crf_struct.mutant_m_sitVisResponse_mean, 2, 'omitnan');
mutant_m_se = std(crf_struct.mutant_m_sitVisResponse_mean, 0, 2, 'omitnan')/sqrt(size(crf_struct.mutant_m_sitVisResponse_mean, 2));

control_f_mean = mean(crf_struct.control_f_sitVisResponse_mean, 2, 'omitnan');
control_f_se = std(crf_struct.control_f_sitVisResponse_mean, 0, 2, 'omitnan')/sqrt(size(crf_struct.control_f_sitVisResponse_mean, 2));

het_f_mean = mean(crf_struct.het_f_sitVisResponse_mean, 2, 'omitnan');
het_f_se = std(crf_struct.het_f_sitVisResponse_mean, 0, 2, 'omitnan')/sqrt(size(crf_struct.het_f_sitVisResponse_mean, 2));

contrasts = ([1, 2, 5, 10, 20, 50, 100]);
test = figure(); hold on;
errorbar(contrasts, control_m_mean, control_m_se, 'k'); plot(contrasts, control_m_mean, 'k.', 'MarkerSize', 12)
errorbar(contrasts, mutant_m_mean, mutant_m_se, 'r'); plot(contrasts, mutant_m_mean, 'r.', 'MarkerSize', 12)
errorbar(contrasts, control_f_mean, control_f_se, 'cyan'); plot(contrasts, control_f_mean, 'cyan.', 'MarkerSize', 12)
errorbar(contrasts, het_f_mean, het_f_se, 'g'); plot(contrasts, het_f_mean, 'g.', 'MarkerSize', 12)
ylim([-0.003 0.012]); xlabel('Contrast (%)'); ylabel('mean V1 response'); 
title('CRF - Not Loco - Males'); 
exportgraphics(test, [fig_directory, 'crf_cdkl5_m_sitting_1s_V1_newnew.pdf'],'ContentType','vector')
