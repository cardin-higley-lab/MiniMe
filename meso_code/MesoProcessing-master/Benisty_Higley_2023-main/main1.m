%% main script for extracting diffusion embedding of correlation traces and modeling behavior
addpath(genpath('D:\meso_demo\MesoProcessing-master\Benisty_Higley_2023-main'))

% processes whole cortex.. really only need one hemisphere


% simulation parameters
winhopSec=0.1; %hop size = # of samples b/w each successive fast fourier transform window
Kfolds = 10; %number of cross validations
winLen = 3; %window size for instantaneous corr mats
J = 20; % how many components 

% load imaging and behavior data
load('K:\rachel\CDKL5\processed_10192023_newWheel\CDKL5\ro138-360\360_12162022\final_timestamps.mat')
load('K:\rachel\CDKL5\processed_10192023_newWheel\CDKL5\ro138-360\360_12162022\Ca_traces_spt_patch9_lowface10p_LSSC.mat')
load('K:\rachel\CDKL5\processed_10192023_newWheel\CDKL5\ro138-360\360_12162022\360_12152022-12162022130722_proc.mat')

nonnancols = find(all(~isnan(parcels_time_trace), 1));
nonnancols(nonnancols>length(spike2_data.greenOnTimestamps)) = [];

imaging_data1.time = spike2_data.greenOnTimestamps(nonnancols);
imaging_data1.data = parcels_time_trace(:,nonnancols);
imaging_data1.fsample = 10;                     
imaging_start_time = imaging_data1.time(1);

behavior_start_ind = findClosestDouble(spike2_data.pupilFrameOnTimestamps, imaging_start_time);
behavior_data.face_time = spike2_data.pupilFrameOnTimestamps(behavior_start_ind:end);
pc1 = proc.motSVD{1, 2}(:,1);
behavior_data.facemap = pc1(behavior_start_ind:end);
behavior_data.pupil = proc.pupil.area(behavior_start_ind:end);

wheel_start_ind = findClosestDouble(spike2_data.analog_signal_time_vect, imaging_start_time);
behavior_data.wheeltime = spike2_data.analog_signal_time_vect(wheel_start_ind:end);
behavior_data.wheel_speed = spike2_data.wheelSpeed(wheel_start_ind:end);



%%

% extract instantaneous correlation matrices (output is 30 less bc that is the size of the window)
[tC, ~] = sliding_window_corr_reg(imaging_data1, winLen, winhopSec);
[sliding_mean, t_win] = sliding_window_mean(imaging_data1, winLen, winhopSec);
t_win=round(t_win);

minLength = min([length(behavior_data.face_time), length(behavior_data.facemap)]);
face_resampled = interp1(behavior_data.face_time(1:minLength), double(behavior_data.facemap(1:minLength)), imaging_data1.time(t_win));
pupil_resampled = interp1(behavior_data.face_time(1:minLength), double(behavior_data.pupil(1:minLength)), imaging_data1.time(t_win));
wheel_resampled = interp1(behavior_data.wheeltime, double(behavior_data.wheel_speed), imaging_data1.time(t_win));


behavior_traces{1} = pupil_resampled;
behavior_traces{2} = face_resampled;
behavior_traces{3} = wheel_resampled;


%%
tic
[recon_signals_corr, R2_corr] = predict_behavior_from_corr(t_win, Kfolds, tC, imaging_data1.time, behavior_traces, J);
toc

tic
[recon_signals_activity, R2_activity] = predict_behavior_from_activity(t_win, Kfolds, sliding_mean, imaging_data1.time, behavior_traces, J);
toc   

ttls = {'pupil','facemap','wheel'};
figure;
for i = 1:3
subplot(3,1,i); plot(imaging_data1.time(t_win), behavior_traces{i});
hold all;
plot(imaging_data1.time(t_win), recon_signals_activity{i});
plot(imaging_data1.time(t_win), recon_signals_corr{i});
ylabel(ttls{i});axis tight;
end
xlabel('Time [sec]');
legend('Behavior','\Phi_a', '\Phi_c');






mask = tril(ones(size(tC,1),size(tC,2)), -1); % make a mask 
B = bsxfun(@times,tC,mask); 


temp = reshape(B, size(tC,1)^2, []);
b = temp(any(temp,2),:);
 

% sort based on std and plot
tc_std = std(b, [], 2)/sqrt(size(b, 1));
[~, tc_order_std] = sort(tc_std, 'ascend');

new_temp = b(tc_order_std, :);

figure();imagesc(new_temp)
hold on; colormap('jet'); clim([0 1]); title('sorted by std')


% sort by mean and plot
tc_means = mean(b, 2);
[~, tc_order_means] = sort(tc_means, 'descend');

new_temp = b(tc_order_means, :);

figure();imagesc(new_temp)
hold on; colormap('jet'); clim([0 1]); title('sorted by mean')



