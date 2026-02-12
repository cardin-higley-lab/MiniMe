function cluster_network_analysis(id_date_time)
%% main script for extracting diffusion embedding of correlation traces and modeling behavior

addpath(genpath('../MesoProcessing-master'))

run('../defineIODirs.m'); % input and output directories

inputDir = fullfile(fixedOutputDir, id_date_time);
outputDir = fullfile('/vast/palmer/scratch/cardin/rlo4/networkOutput/', id_date_time);

%inputDir ='K:\rachel\CDKL5\processed_01242024\CDKL5\ro138-360\360_12162022'; outputDir = inputDir;

mkNewDir(outputDir)

% simulation parameters
winhopSec=0.1; %hop size = # of samples b/w each successive fast fourier transform window
Kfolds = 10; %number of cross validations
winLen = 3; %window size for instantaneous corr mats
J = 20; % how many components 

%% 

load(fullfile(inputDir, 'facecam_motion_energy.mat'), 'vid_energy');
load(fullfile(inputDir, 'final_timestamps.mat'),'spike2_data');
load(fullfile(inputDir, 'Ca_traces_spt_patch9_lowface10p_LSSC.mat'), 'parcels_time_trace');

if isfield(spike2_data, 'blueOnTimestamps')
    mesoTimestampsOn = spike2_data.blueOnTimestamps;
else
    mesoTimestampsOn = spike2_data.greenOnTimestamps;
end


half_parcels = parcels_time_trace(1:size(parcels_time_trace, 1)/2, :);
nonnancols = find(all(~isnan(half_parcels), 1));
nonnancols(nonnancols>length(mesoTimestampsOn)) = [];

imaging_data1.time = mesoTimestampsOn(nonnancols);
imaging_data1.data = half_parcels(:,nonnancols);
imaging_data1.fsample = 10;                     
imaging_start_time = imaging_data1.time(1);

behavior_start_ind = findClosestDouble(spike2_data.pupilFrameOnTimestamps, imaging_start_time);
behavior_data.face_time = spike2_data.pupilFrameOnTimestamps(behavior_start_ind:end);
behavior_data.facemap = vid_energy(behavior_start_ind:end);

wheel_start_ind = findClosestDouble(spike2_data.analog_signal_time_vect, imaging_start_time);
behavior_data.wheeltime = spike2_data.analog_signal_time_vect(wheel_start_ind:end);
behavior_data.wheel_speed = spike2_data.wheelSpeed(wheel_start_ind:end);



%%

% extract instantaneous correlation matrices (output is 30 less bc that is the size of the window)
[tC, ~] = sliding_window_corr_reg(imaging_data1, winLen, winhopSec);
save(fullfile(outputDir, 'tCorrs.mat'), 'tC', '-v7.3'); 

[sliding_mean, t_win] = sliding_window_mean(imaging_data1, winLen, winhopSec);
t_win=round(t_win);

minLength = min([length(behavior_data.face_time), length(behavior_data.facemap)]);
face_resampled = interp1(behavior_data.face_time(1:minLength), double(behavior_data.facemap(1:minLength)), imaging_data1.time(t_win));
wheel_resampled = interp1(behavior_data.wheeltime, double(behavior_data.wheel_speed), imaging_data1.time(t_win));


behavior_traces{1} = face_resampled;
behavior_traces{2} = wheel_resampled;

%% calculate diffmap and save

if exist(fullfile(outputDir, 'riemmaniam_mean.mat'), 'file') == 2 % 2 checks for a file
    disp('riemmaniam mean already calculated')
    load(fullfile(outputDir, 'riemmaniam_mean.mat'),  'mX');
else
    % get riemmannian projections for instantaneous correlations
    disp('calculating riemmaniam mean')
    mRiemannianMean = RiemannianMean(tC);
    mX = proj_R1(mRiemannianMean^(-1/2), tC);
    save(fullfile(outputDir, 'riemmaniam_mean.mat'), 'mX', '-v7.3');
end


if exist(fullfile(outputDir, 'diffmap.mat'), 'file') == 2 % 2 checks for a file
    disp('diffmap already calculated')
    load(fullfile(outputDir, 'diffmap.mat'), 'diffmap');
else
    disp('calculating initAll')
    par.knn = max(round(size(mX,2)*0.01), 20);
    [ initAll, ~ ] = CalcInitAff2D(mX, par );
    
    disp('calculating diffmap')
    dParams.maxInd = min(size(initAll,1), 1+J);
    [~, ~, Psi] = calcDiffusionMap(initAll, dParams);
    diffmap = Psi(:,2:end).';
    save(fullfile(outputDir, 'diffmap.mat'), 'diffmap', '-v7.3');
end





%%
tic
[recon_signals_corr, R2_corr] = predict_behavior_from_corr(t_win, Kfolds, tC, diffmap, imaging_data1.time, behavior_traces, J);
toc
save(fullfile(outputDir, 'recon_signals_corr.mat'), 'recon_signals_corr', 'R2_corr', '-v7.3'); 

%load(fullfile(outputDir, 'recon_signals_corr.mat'),'recon_signals_corr', 'R2_corr');

tic
[recon_signals_activity, R2_activity] = predict_behavior_from_activity(t_win, Kfolds, sliding_mean, imaging_data1.time, behavior_traces, J);
toc   
save(fullfile(outputDir, 'recon_signals_activity.mat'), 'recon_signals_activity', 'R2_activity', '-v7.3'); 


%% full model

[res_xraw_phi_ridge, res_xrawsh_phi_ridge, res_xraw_phish_ridge, predicted_xraw_phi_ridge, predicted_xrawsh_phi_ridge, predicted_xraw_phish_ridge]= train_test_modeling_behavior(sliding_mean, imaging_data1.time, behavior_traces, diffmap, t_win, J, Kfolds);
save(fullfile(outputDir, 'combinedModel.mat'), 'res_xraw_phi_ridge', 'res_xrawsh_phi_ridge', 'res_xraw_phish_ridge', 'predicted_xraw_phi_ridge',...
    'predicted_xrawsh_phi_ridge', 'predicted_xraw_phish_ridge', '-v7.3'); 


%% figs 

ttls = {'facemap','wheel'};
model_res = figure;
for i = 1:length(ttls)
    subplot(length(ttls),1,i); plot(imaging_data1.time(t_win), behavior_traces{i});
    hold all;
    plot(imaging_data1.time(t_win), recon_signals_activity{i}); xlim([200 1200])
    plot(imaging_data1.time(t_win), recon_signals_corr{i}); xlim([200 1200])
    ylabel(ttls{i});axis tight;
end
xlabel('Time [sec]');
legend('Behavior','Phi a', 'Phi c');

saveas(model_res, fullfile(outputDir, 'model_res')); 
%exportgraphics(model_res,'C:\Users\CardinLab\Desktop\360_12162022\model_res.pdf','ContentType','vector')



%% more figs

mask = tril(ones(size(tC,1),size(tC,2)), -1); % make a mask 
B = bsxfun(@times,tC,mask); 


temp = reshape(B, size(tC,1)^2, []);
b = temp(any(temp,2),:);
 

% sort based on std and plot
tc_std = std(b, [], 2)/sqrt(size(b, 1));
[~, tc_order_std] = sort(tc_std, 'ascend');

new_temp = b(tc_order_std, :);

std_waterfall = figure();imagesc(abs(new_temp))
hold on; colormap('jet'); clim([0 1]); title('sorted by std')

saveas(std_waterfall, fullfile(outputDir, 'std_waterfall')); 




% sort by mean and plot
tc_means = mean(b, 2);
[~, tc_order_means] = sort(tc_means, 'descend');

new_temp = b(tc_order_means, :);

mean_waterfall = figure();imagesc(new_temp)
hold on; colormap('jet'); clim([0 1]); title('sorted by mean')

saveas(mean_waterfall, fullfile(outputDir, 'mean_waterfall')); 



% mean_waterfall = figure();imagesc(new_temp)
% hold on; colormap('jet'); clim([0 1]); title('sorted by mean'); xlim([0 2000]); colorbar();
% exportgraphics(mean_waterfall,'C:\Users\CardinLab\Desktop\360_12162022\mean_waterfall.pdf','ContentType','vector')



end
















