%% CRF analysis

clear vars

dir_path = 'K:\rachel\CRF\processed_09062023\';
masterDir = dir(dir_path);


fs = 10; preSeconds = 1; postSeconds = 5;
t = -preSeconds:(preSeconds+postSeconds)/(preSeconds*fs+postSeconds*fs):postSeconds;


for animals = 1:length(masterDir)

    animal = masterDir(animals).name;
    if length(animal) < 3
        continue;
    end

    subfolders = dir(fullfile(dir_path, animal));
    for sessions = 1:length(subfolders)

        session = subfolders(sessions).name;
        if length(session) < 3
            continue;
        end

        folder = fullfile(dir_path, animal, session);
        disp(['working on: ' folder])

        load(fullfile(folder, 'final_timestamps.mat'), 'spike2_data');
        load(fullfile(folder, 'behavior_timestamps.mat'), 'behavior');
        load(fullfile(folder, 'Ca_traces_spt_patch9_lowface10p_Allen.mat'), 'parcels_time_trace');

        vis_folder = fullfile(folder, 'vis_data');
        mkNewDir(vis_folder);
        
        if ~exist(fullfile(folder, 'stim_order.mat'), "file") == 2
            xcelList = (dir(fullfile(folder, '*.csv')));
            excelFile = fullfile(folder,xcelList.name);
            crf = readtable(excelFile);
            contrasts = crf.Var1 ;
            contrast_vals = unique(contrasts);
            save(fullfile(folder, 'stim_order.mat'), 'contrasts');
        else
            load(fullfile(folder, 'stim_order.mat'), 'contrasts')
            contrast_vals = unique(contrasts);
        end

        if length(contrasts) ~= length(spike2_data.diodeOffTimestamps)
            warning('length of vis stims doesnt match diode timestamps')
        end

        [vis_responses, loco] = calc_crf_responses(contrasts, parcels_time_trace, spike2_data, behavior);

        % crf_plot = figure(); hold on;
        % v1_response = mean(vis_responses(2,:,:), 3, 'omitnan');
        % v1_response_se = std(vis_responses(2,:,:), 0, 3, 'omitnan')/sqrt(size(vis_responses, 3));
        % shadedErrorBar(contrast_vals, v1_response, v1_response_se, 'lineprops', {'k', 'markerfacecolor','k'});
        % plot((contrast_vals), v1_response, 'ko')
        % xlabel('contrast value'); ylabel('mean v1 response'); title('Contrast response function')
        % exportgraphics(crf_plot, fullfile(vis_folder, 'crf_plot.pdf'),'ContentType','vector')

        % get loco and sit-specific responses
        loco_responses = nan(size(vis_responses));
        loco_responses(:, loco==1) = vis_responses(:, loco==1);
        mean_loco_responses = mean(loco_responses, 3, 'omitnan');
        se_loco_responses = std(loco_responses, 0, 3, 'omitnan')/sqrt(size(loco_responses, 3));

        % crf_plot = figure(); hold on;
        % shadedErrorBar(contrast_vals, mean_loco_responses(1,:), se_loco_responses(1,:), 'lineprops', {'k', 'markerfacecolor','k'});
        % plot((contrast_vals), mean_loco_responses(1,:), 'ko')
        % xlabel('contrast value'); ylabel('mean v1 response'); title('Contrast response function - loco')
        % exportgraphics(crf_plot, fullfile(vis_folder, 'crf_plot_loco.pdf'),'ContentType','vector')

        sit_responses = nan(size(vis_responses));
        sit_responses(:, loco==0) = vis_responses(:, loco==0);
        mean_sit_responses = mean(sit_responses, 3, 'omitnan');
        se_sit_responses = std(sit_responses, 0, 3, 'omitnan')/sqrt(size(sit_responses, 3));
    
        % crf_plot = figure(); hold on;
        % shadedErrorBar(contrast_vals, mean_sit_responses(1,:), se_sit_responses(1,:), 'lineprops', {'k', 'markerfacecolor','k'});
        % plot((contrast_vals), mean_sit_responses(1,:), 'ko')
        % xlabel('contrast value'); ylabel('mean v1 response'); title('Contrast response function - sitting')
        % exportgraphics(crf_plot, fullfile(vis_folder, 'crf_plot_sit.pdf'),'ContentType','vector')

        % %get timestamps for all contrasts
        % v1_100p = spike2_data.diodeOnTimestamps(find(contrasts==5));
        % 
        % [alignedData_loco, preData_loco] = averageMovieToEvent(parcels_time_trace, spike2_data.greenOnTimestamps, v1_100p(loco(7,:)==1), preSeconds*fs, postSeconds*fs);
        % normed_aligned_loco = alignedData_loco - mean(preData_loco, 3);
        % v1_loco_avg = squeeze(mean(normed_aligned_loco(:,1,:), 1, 'omitnan'));
        % 
        % [alignedData_sit, preData_sit] = averageMovieToEvent(parcels_time_trace, spike2_data.greenOnTimestamps, v1_100p(loco(7,:)==0), preSeconds*fs, postSeconds*fs);
        % normed_aligned_sit = alignedData_sit - mean(preData_sit, 3);
        % v1_sit_avg = squeeze(mean(normed_aligned_sit(:,1,:), 1, 'omitnan'));

        % figure(); hold on
        % subplot(2,1,1); plot(t(2:end), v1_sit_avg); title('wheel off'); xline(0); ylim([min([v1_loco_avg; v1_sit_avg]), max([v1_loco_avg; v1_sit_avg])])
        % subplot(2,1,2); plot(t(2:end), v1_loco_avg); title('wheel on'); xline(0); ylim([min([v1_loco_avg; v1_sit_avg]), max([v1_loco_avg; v1_sit_avg])])
        % sgtitle(session,'Interpreter',"none")

        
        [alignedData, preData] = averageMovieToEvent(parcels_time_trace, spike2_data.greenOnTimestamps, spike2_data.diodeOnTimestamps, preSeconds*fs, postSeconds*fs);
        normed_aligned = alignedData./mean(preData, 3);
        v1_aligned_1p = normed_aligned((contrasts==1),1,:);
        v1_aligned_2p = normed_aligned((contrasts==2),1,:);
        v1_aligned_5p = normed_aligned((contrasts==5),1,:);
        v1_aligned_10p = normed_aligned((contrasts==10),1,:);
        v1_aligned_20p = normed_aligned((contrasts==20),1,:);
        v1_aligned_50p = normed_aligned((contrasts==50),1,:);
        v1_aligned_100p = normed_aligned((contrasts==100),1,:);


        %save(fullfile(vis_folder, 'crf_data'), "vis_responses", "loco", "loco_responses", "sit_responses", ...
            %"mean_loco_responses", "se_loco_responses", "mean_sit_responses", "se_sit_responses", ...
            %"normed_aligned_loco", "normed_aligned_sit")

        save(fullfile(vis_folder, 'crf_trace_data'), "v1_aligned_1p", "v1_aligned_2p", "v1_aligned_5p",  ...
            "v1_aligned_10p", "v1_aligned_20p", "v1_aligned_50p", "v1_aligned_100p", "loco", "fs", "preSeconds", "postSeconds")

    end


end