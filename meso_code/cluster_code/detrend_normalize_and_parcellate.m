function detrend_normalize_and_parcellate(id_date_time, norm_method)

fprintf([id_date_time ' started processing at %s\n'], datestr(now,'HH:MM:SS.FFF'))


%id_date_time = 'temp\360_12162022';

addpath(genpath('../MesoProcessing-master'))
load('parcells_updated121519.mat','parcells_new');
load('brain_mask.mat','brain_mask');
inds2keep = find(brain_mask);

run('../defineIODirs.m'); % input and output directories

outputDir = fullfile(fixedOutputDir,id_date_time);
load(fullfile(outputDir, 'tform.mat'),'tform','R','C');
load(fullfile(outputDir, 'raw_mats.mat'),'blue_raw');
load(fullfile(outputDir, 'final_timestamps.mat'),'spike2_data');
load(fullfile(outputDir, 'facecam_motion_energy.mat'), 'vid_energy')

outputDir = fullfile(fixedOutputDir,id_date_time);
findslash = strfind(id_date_time,'/');
animal = id_date_time(findslash(1)+1:findslash(2)-1);
session = id_date_time(findslash(2)+1:end);

 %% reprocess states 

        states = newStates(spike2_data, vid_energy, 2, 3, 10);
        save(fullfile(outputDir, 'state_timestamps.mat'), 'states', '-v7.3');

        statesplot = states_fig(states, spike2_data, vid_energy, animal, session);
        saveas(statesplot, fullfile(outputDir, 'statesplot'));

        behavior = newStates(spike2_data, vid_energy, 2, 0, 0);
        save(fullfile(outputDir, 'behavior_timestamps.mat'), 'behavior', '-v7.3');


%% Alignment 
tic

disp('Aligning to template...');

% transform
blue_aligned = transform_frames(blue_raw, tform, R, C);
clear blue_raw

% save clip of output
filename = fullfile(outputDir, 'alignedOutput.gif');
meso_clip(blue_aligned, filename, 1400:1500, 'smoothed_no')

disp(['Aligning to atlas: ' num2str(toc)])

%% trim and reshape
tic

blue_aligned = single(reshape(blue_aligned, 256, 256, []));
dff_blue = cat(3, nan(size(blue_aligned(:, :, 1:delay_length))), blue_aligned(:, :, delay_length+1:end)); % converts beginning of session to nans
dff_blue = reshape(dff_blue, R*C, size(dff_blue, 3));
output_sig = nan(size(dff_blue));
output_sig(inds2keep, :) = dff_blue(inds2keep, :);

% save full output
% save(fullfile(outputDir, 'full_output.mat'),'output_sig', '-v7.3');


%% normalization

switch norm_method
    case 'none'% do nothing
        normed_sig = output_sig;

    case 'zscored' %redoes zscore and then does dff
        [normed_sig, ~, ~] = zscore_pixelwise(output_sig, delay_length);

    case 'dff' %finds df/f using mean of bottom 10% as f0
        [normed_sig] = dff_bottom10(output_sig);

    case 'lowface10p'
        if isfield(spike2_data, 'blueOnTimestamps')
            mesoTimestampsOn = spike2_data.blueOnTimestamps; 
        else
            mesoTimestampsOn = spike2_data.greenOnTimestamps; 
        end
        [normed_sig, F0, ~] = lowfaceNorm_bottom10(states, mesoTimestampsOn, output_sig);
        f0_fig = figure();
        imagesc(reshape(F0, 256, 256));
        saveas(f0_fig, fullfile(outputDir, 'f0.png'))

    otherwise
        error('Illeagal Normalization method');
end
    
    normalized = fullfile(outputDir, ['normed_' norm_method '_blue_uv.mat']);
    save(normalized, 'normed_sig', '-v7.3'); % save time trace both w and without mean removal

    disp('Normalization complete')

%% save short clip for confirming output

filename = fullfile(outputDir, 'confirmOutput.gif');
meso_clip(normed_sig, filename, 1400:1500, 'smoothed_no')


%% save short clip for confirming output (smoothed)

filename = fullfile(outputDir, 'confirmOutput_sm.gif');
meso_clip(normed_sig, filename, 1400:1500, 'smoothed_yes')


%% Parcellation: 'allen'

    brain_mask_new = brain_mask; 
    brain_mask_new(:, 119:138) = false; 
    brain_mask_new = reshape(brain_mask_new, 256*256, []);
    normed_sig(~brain_mask_new, :) = NaN;
    parcels_time_trace = parcellate_allen(parcells_new, normed_sig);

    allenfile = fullfile(outputDir,  ['Ca_traces_spt_patch9_' norm_method '_Allen.mat']); %find me: _Allen_dffnorm.mat
    save(allenfile, 'parcels_time_trace'); % save time trace both w and without mean removal

    disp('Allen parcellation complete')

%% Parcellation: 'LSSC'

    % remove the mid-line
    [R,C] = size(brain_mask); 
    brain_mask(:, 121:136) = false; %previously: 125:131
    brain_mask_left = brain_mask;
    brain_mask_right = brain_mask;
    brain_mask_left(:,129:end) = false;
    brain_mask_right(:,1:128) = false;
   
    % we parcellate left and right separately
    cfg.outputfilepath = outputDir;
    cfg.preProcess=false;
    cfg.N_TRIALS=10;
    cfg.n_clust = 5000 ;
    cfg.makePlots = false;
    cfg.thrcluster=0.99; % change this to influence number of parcels; can give a vector of different threshold (e.g. [0.9:0.01:0.99])
    cfg.NROWS = R;
    cfg.NCOLS = C;
    cfg.isoverlap = false;
    
    allregionspix_left = find(brain_mask_left);
    dFoF_masked_left = normed_sig(allregionspix_left, :); 
    st=1e3;en=size(dFoF_masked_left,2)-1e3;
    cfg.title_str ='left_LSSC';
    warning('off', 'MATLAB:eigs:IgnoredOptionIssym');
    runROI_meso_nlm(cfg, dFoF_masked_left(:,st:en), allregionspix_left, brain_mask_left);
    clear dFoF_masked_left brain_mask_left allregionspix_left

    allregionspix_right = find(brain_mask_right);
    dFoF_masked_right = normed_sig(allregionspix_right, :); 
    st=1e3;en=size(dFoF_masked_right,2)-1e3;
    cfg.title_str ='right_LSSC';
    runROI_meso_nlm(cfg, dFoF_masked_right(:,st:en), allregionspix_right, brain_mask_right);
    clear dFoF_masked_right brain_mask_right allregionspix_right
    
    leftfile = dir(fullfile(outputDir, 'left_LSSC*.mat'));
    parleft = load(fullfile(leftfile.folder, leftfile.name));
    rightfile = dir(fullfile(outputDir, 'right_LSSC*.mat'));
    parright = load(fullfile(rightfile.folder, rightfile.name));
    parcelsAll = [parleft.ROI_list parright.ROI_list];
   
    
    parmap = zeros(R,C);
    parcels_time_trace = zeros(length(parcells_new.names), size(normed_sig,2), 'single');  
    
    for par_i = 1:length(parcelsAll)
        roiinds=parcelsAll(par_i).pixel_list;
        parmap(roiinds) = par_i;
        
        pardata = normed_sig(roiinds, :);        
        parcels_time_trace(par_i,:) = mean(pardata, 'omitnan');
         
    end
    
    galfile = fullfile(outputDir,  ['Ca_traces_spt_patch9_' norm_method '_LSSC.mat']); %find me: _Allen_dffnorm.mat
    save(galfile, 'parcels_time_trace', 'parmap'); % 'parcels_time_trace_GSR', 

    delete(fullfile(leftfile.folder, leftfile.name))
    delete(fullfile(rightfile.folder, rightfile.name))
    close all

    disp(['deleting unnecessary files: ' num2str(toc) 's'])

    disp('LSSC parcellation complete')


%% Parcellation: Grids

load("grid_obj.mat", 'grid_obj')
num_grid_parcels = size(grid_obj.indicators, 2);

parcels_time_trace = zeros(num_grid_parcels, size(normed_sig,2), 'single');
for par_i = 1:num_grid_parcels

    roiinds = grid_obj.indicators(:,par_i) == 1;
    pardata = normed_sig(roiinds, :);
    detect = isoutlier(pardata);
    keep_data = nan(size(pardata));
    keep_data(detect==0) = pardata(detect==0);
    parcels_time_trace(par_i,:) = mean(keep_data, 'omitnan');

end

gridfile = fullfile(outputDir,  ['Ca_traces_spt_patch9_' norm_method '_Grid4.mat']);
save(gridfile, 'parcels_time_trace');  

disp('grids parcellation complete')


%% sanity figs

% define tings and load things needed for making figs

load(allenfile, 'parcels_time_trace');
load(fullfile(outputDir, 'facecam_motion_energy.mat'), 'vid_energy');

if isfield(spike2_data, 'blueOnTimestamps')
    mesoTimestampsOn = spike2_data.blueOnTimestamps;
else
    mesoTimestampsOn = spike2_data.greenOnTimestamps;
end

wheel_fs = 5000; imaging_fs = 10; preSeconds = 2; postSeconds = 6;
locoOn = states.locoOn-3;
whiskOn = states.faceHighSitOn;


% wheel and meso alignment check
sanityCheck = figure();
ax(1) = subplot(211);
plot(spike2_data.analog_signal_time_vect, spike2_data.wheelSpeed);
ylabel('wheel speed'); title(['alignment check - ' id_date_time], 'Interpreter', 'none');
ax(2) = subplot(212);
maxLength = min(length(mesoTimestampsOn), size(parcels_time_trace, 2));
plot(mesoTimestampsOn(1:maxLength), mean(parcels_time_trace(:, 1:maxLength),1, 'omitnan'));
ylabel('mean(all parcels)'); xlabel('LEDonTimestamps (seconds)'); linkaxes(ax,'x');
saveas(sanityCheck, fullfile(outputDir, 'alignmentCheck'));




% getting data for state-dep activity
if size(spike2_data.wheelSpeed, 1) > 1; spike2_data.wheelSpeed = spike2_data.wheelSpeed'; end
[alignedData, ~] = averageMovieToEvent(spike2_data.wheelSpeed, spike2_data.analog_signal_time_vect, locoOn, preSeconds*wheel_fs, postSeconds*wheel_fs);
aligned_zscored = zscore(squeeze(alignedData), [], 2);
baseline_means = mean(aligned_zscored(:,1:preSeconds*wheel_fs), 2);
aligned_zscored_baselinesub = aligned_zscored - baseline_means;
avgWheelOnset = squeeze(mean(aligned_zscored_baselinesub,1, 'omitnan'));

[alignedData, ~] = averageMovieToEvent(vid_energy', spike2_data.pupilFrameOnTimestamps, whiskOn, preSeconds*imaging_fs, postSeconds*imaging_fs);
aligned_zscored = zscore(squeeze(alignedData), [], 2);
if size(aligned_zscored, 2)>1
    baseline_means = mean(aligned_zscored(:,1:preSeconds*imaging_fs), 2);
    aligned_zscored_baselinesub = aligned_zscored - baseline_means;
    avgWhiskOnset = squeeze(mean(aligned_zscored_baselinesub,1, 'omitnan'));
else
    baseline_means = mean(aligned_zscored(1:preSeconds*imaging_fs, :), 1);
    aligned_zscored_baselinesub = aligned_zscored' - baseline_means;
    avgWhiskOnset = squeeze(mean(aligned_zscored_baselinesub,1, 'omitnan'));
end

[alignedData, ~] = averageMovieToEvent(parcels_time_trace(1,:), mesoTimestampsOn, locoOn, preSeconds*imaging_fs, postSeconds*imaging_fs);
aligned_zscored = zscore(squeeze(alignedData), [], 2);
baseline_means = mean(aligned_zscored(:,1:preSeconds*imaging_fs), 2);
aligned_zscored_baselinesub = aligned_zscored - baseline_means;
avgLocoActivity = squeeze(mean(aligned_zscored_baselinesub,1, 'omitnan'));

[alignedData, ~] = averageMovieToEvent(parcels_time_trace(1,:), mesoTimestampsOn, whiskOn, preSeconds*imaging_fs, postSeconds*imaging_fs);
aligned_zscored = zscore(squeeze(alignedData), [], 2);
if size(aligned_zscored, 2)>1
    baseline_means = mean(aligned_zscored(:,1:preSeconds*imaging_fs), 2);
    aligned_zscored_baselinesub = aligned_zscored - baseline_means;
    avgWhiskActivity = squeeze(mean(aligned_zscored_baselinesub,1, 'omitnan'));
else
    baseline_means = mean(aligned_zscored(1:preSeconds*imaging_fs, :), 1);
    aligned_zscored_baselinesub = aligned_zscored' - baseline_means;
    avgWhiskActivity = squeeze(mean(aligned_zscored_baselinesub,1, 'omitnan'));
end


% state-dep activity plot
state_activity_plot = figure(); hold on
ax(1) = subplot(2, 2, 1); plot(avgWheelOnset); title('avg loco Onset'); xline(10000); ylim([-1 4])
ax(2) = subplot(2, 2, 2); plot(avgWhiskOnset); title('avg whisk Onset'); xline(20); ylim([-1 4])
ax(3) = subplot(2, 2, 3); plot(avgLocoActivity); title('v1 activity at loco onset'); xline(20); ylim([-1 4])
ax(4) = subplot(2, 2, 4); plot(avgWhiskActivity); title('v1 activity at whisk onset'); xline(20); ylim([-1 4])
saveas(state_activity_plot, fullfile(outputDir, 'state_check.fig'));



%% done

disp('sanity figures complete')
fprintf([id_date_time ' finished processing at %s\n'], datestr(now,'HH:MM:SS.FFF'))

end