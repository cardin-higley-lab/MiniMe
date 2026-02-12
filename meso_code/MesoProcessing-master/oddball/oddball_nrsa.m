
addpath(genpath('D:\meso_demo\analysis_scripts'))
addpath(genpath('D:\meso_demo\MesoProcessing-master'))


session = 'C:\Users\CardinLab\Desktop\Oddball_testing\59_platformtest\59_oddball_piezo_10122022';
mkNewDir(fullfile(session, 'oddball'))
cd(fullfile(session, 'oddball'))

%% load neural signal 

load(fullfile(session, 'Ca_traces_spt_patch9_Allen.mat'))
zData = (parcels_time_trace-mean(parcels_time_trace,2, 'omitnan'))./std(parcels_time_trace,[],2, 'omitnan');
v1_data = mean([zData(1,:); zData(2,:)]);


%% load oddball
load(fullfile(session, 'oddballStimuliData.mat'));

%% load and process spike2 
load(fullfile(session, 'final_timestamps.mat'));
mergeOn = spike2_data.diodeOnTimestamps;
mergeOff = spike2_data.diodeOffTimestamps;
time_diffs = mergeOff-mergeOn;
false_times = find(time_diffs> 0.6);
if ~isempty(false_times) 
    mergeOn(false_times) = [];
    mergeOff(false_times) = [];
end

if length(mergeOn) ~= length(oddballStimuliData.masterSequence)
    warning('unexpected number of diode pulses')
    disp(['number of diode pulses: ' num2str(length(mergeOn))])
end

%% load facemap data
% facemapfile = ls(fullfile(session, '\*proc.mat'));
% load(fullfile(session, facemapfile), 'proc')
% face_proc = proc;
% 
% pc = face_proc.motSVD{1, 2}(:,1);
% if skewness(pc) < 0.5; pc = -pc; end

%% define and trim data

%maxlength = min(min(length(spike2_data.pupilFrameOnTimestamps), length(pc)), min(size(zData,2), length(spike2_data.blueOnTimestamps)));
maxlength = min(length(spike2_data.pupilFrameOnTimestamps), min(size(zData,2), length(spike2_data.blueOnTimestamps)));
pupil_time = spike2_data.pupilFrameOnTimestamps(1:maxlength);
%face_Norm = normalize(pc(1:maxlength));
zData = zData(:, 1:maxlength);
mesoTimestamps = spike2_data.blueOnTimestamps(1:maxlength);


%% define face states

% load(fullfile(session, 'state_timestamps.mat'))
% sitOn = mesoTimestamps(1);
% sitOff = mesoTimestamps(end);
% 
% 
% % gets time points within sitting bouts when not grooming
% [notGroomOn, notGroomOff] = timestampsNot(states.groomHighOn, states.groomHighOff, 'StartTime', sitOn, 'EndTime', sitOff);
% [sitOn_final, sitOff_final] = timestampsAnd(sitOn, sitOff, notGroomOn, notGroomOff);
% 
% % gets time points within sitting bouts 
% pupilTime_Idx=cell(1,length(sitOn_final));
% for st=1:length(sitOn_final)
%     pupilTime_Idx{st}=find(pupil_time > sitOn_final(st) & pupil_time < sitOff_final(st));
% end
% pupilTime=cell2mat(pupilTime_Idx');
% face=face_Norm(pupilTime);
% 
% % define quantiles 
% zThresh_bottom1 = quantile(face,0.75);
% zThresh_bottom2 = quantile(face,0.5);
% zThresh_bottom3 = quantile(face,0.25);


%% get oddball experiment parameters

manyControlsLength = oddballStimuliData.numPresentations;

epochLengths = [oddballStimuliData.numPresentations oddballStimuliData.preSequenceRedundantStimuli oddballStimuliData.numPresentations oddballStimuliData.preSequenceRedundantStimuli oddballStimuliData.numPresentations];
epochMasks = false(sum(epochLengths),length(epochLengths));
currStart = 0;
for ii = 1:length(epochLengths)
  epochMasks(currStart+1:currStart+epochLengths(ii),ii) = true;
  currStart = currStart + epochLengths(ii);
end

% let's define the first stim 
stimA = mode(oddballStimuliData.masterSequence(epochMasks(:,5))); % A = oddball first

controlstimATimestampsOn = mergeOn(epochMasks(:,1) & oddballStimuliData.masterSequence==stimA);
controlstimATimestampsOff = mergeOff(epochMasks(:,1) & oddballStimuliData.masterSequence==stimA);
stimA_dev_TimestampsOn = mergeOn(epochMasks(:,3) & oddballStimuliData.masterSequence==stimA);
stimA_dev_TimestampsOff = mergeOff(epochMasks(:,3) & oddballStimuliData.masterSequence==stimA);
stimA_preBlock3TimestampsOn = mergeOn(epochMasks(:,4));
stimA_preBlock3TimestampsOff = mergeOff(epochMasks(:,4));
stimA_red_TimestampsOn = mergeOn(epochMasks(:,5) & oddballStimuliData.masterSequence==stimA);
stimA_red_TimestampsOff = mergeOff(epochMasks(:,5) & oddballStimuliData.masterSequence==stimA);

% ok now the other stim
stimB = mode(oddballStimuliData.masterSequence(epochMasks(:,3))); % B = oddball second

controlstimBTimestampsOn = mergeOn(epochMasks(:,1) & oddballStimuliData.masterSequence==stimB);
controlstimBTimestampsOff = mergeOff(epochMasks(:,1) & oddballStimuliData.masterSequence==stimB);
stimB_preBlock2TimestampsOn = mergeOn(epochMasks(:,2)); 
stimB_preBlock2TimestampsOff = mergeOff(epochMasks(:,2));
stimB_red_TimestampsOn = mergeOn(epochMasks(:,3) & oddballStimuliData.masterSequence==stimB);
stimB_red_TimestampsOff = mergeOff(epochMasks(:,3) & oddballStimuliData.masterSequence==stimB);
stimB_dev_TimestampsOn = mergeOn(epochMasks(:,5) & oddballStimuliData.masterSequence==stimB);
stimB_dev_TimestampsOff = mergeOff(epochMasks(:,5) & oddballStimuliData.masterSequence==stimB);




%% LINEGRAPHS WITH ERROR BARS

% looking at both blocks together
preSeconds = 1; postSeconds = 2; fs=10;
deviant_traces = averageMovieToEvent(zData,mesoTimestamps,[stimA_dev_TimestampsOn; stimB_dev_TimestampsOn],preSeconds*fs,postSeconds*fs);
deviant_traces_v1 = squeeze(mean(deviant_traces(:, [1;2], :), 2, 'omitnan'));
deviant_v1_avgTrace = mean(deviant_traces_v1, 1, 'omitnan');
deviant_v1_se = std(deviant_traces_v1, 0, 1, 'omitnan')/sqrt(size(deviant_traces_v1, 1));

redundant_traces = averageMovieToEvent(zData,mesoTimestamps,[stimA_red_TimestampsOn; stimB_red_TimestampsOn],10,20);
redundant_traces_v1 = squeeze(mean(redundant_traces(:, [1;2], :), 2, 'omitnan'));
redundant_v1_avgTrace = mean(redundant_traces_v1, 1, 'omitnan');
redundant_v1_se = std(redundant_traces_v1, 0, 1, 'omitnan')/sqrt(size(redundant_traces_v1, 1));

control_traces = averageMovieToEvent(zData,mesoTimestamps,[mergeOn(1:manyControlsLength)],10,20);
control_traces_v1 = squeeze(mean(control_traces(:, [1;2], :), 2, 'omitnan'));
control_v1_avgTrace = mean(control_traces_v1, 1, 'omitnan');
control_v1_se = std(control_traces_v1, 0, 1, 'omitnan')/sqrt(size(control_traces_v1, 1));

test = figure(); hold on;
t = -preSeconds+preSeconds/fs:(preSeconds+postSeconds)/(preSeconds*fs+postSeconds*fs):postSeconds; % get time vect for plotting
control = shadedErrorBar(t, control_v1_avgTrace, control_v1_se, 'lineprops', {'b', 'markerfacecolor','b'}); 
redundant = shadedErrorBar(t, redundant_v1_avgTrace, redundant_v1_se, 'lineprops', {'k', 'markerfacecolor','k'}); 
deviant = shadedErrorBar(t, deviant_v1_avgTrace, deviant_v1_se, 'lineprops', {'r', 'markerfacecolor','r'}); 
xlabel('time (s)'); ylabel('dff')
xline(0)
legend([control.patch redundant.patch deviant.patch], 'avg control response', 'avg redundant response', 'avg deviant response', 'Location', 'northeast');
exportgraphics(test, 'linegraph_errorbars.pdf','ContentType','vector')



save("linegraph_data.mat", "t", "control_v1_avgTrace", "control_v1_se", "redundant_v1_avgTrace", "redundant_v1_se", "deviant_v1_avgTrace", "deviant_v1_se")


% % block 1
% block1_dev_v1_traces = deviant_traces_v1(1:length(stimA_dev_TimestampsOn), :);
% block1_dev_v1_avgTrace = mean(block1_dev_v1_traces, 1, 'omitnan');
% block1_red_v1_traces = redundant_traces_v1(1:length(stimB_red_TimestampsOn), :);
% block1_red_v1_avgTrace = mean(block1_red_v1_traces, 1, 'omitnan');
% 
% figure(); hold on
% plot(block1_red_v1_avgTrace, 'k')
% plot(block1_dev_v1_avgTrace, 'r')
% legend({'avg redundant response', 'avg deviant response'}, 'Location', 'northwest')
% ylim([-0.4 1])
% 
% % block 2
% block2_dev_v1_traces = deviant_traces_v1(length(stimB_dev_TimestampsOn)+1:end, :);
% block2_dev_v1_avgTrace = mean(block2_dev_v1_traces, 1, 'omitnan');
% block2_red_v1_traces = redundant_traces_v1(length(stimA_red_TimestampsOn)+1:end, :);
% block2_red_v1_avgTrace = mean(block2_red_v1_traces, 1, 'omitnan');
% 
% figure(); hold on
% plot(block2_red_v1_avgTrace, 'k')
% plot(block2_dev_v1_avgTrace, 'r')
% legend({'avg redundant response', 'avg deviant response'}, 'Location', 'northwest')
% ylim([-0.4 1])


%% barplot for responses 

control_responses = mean(control_traces_v1(:,preSeconds*fs:fs*2), 2) - mean(control_traces_v1(:,fs/2+1:fs), 2); %
red_responses = mean(redundant_traces_v1(:,preSeconds*fs:fs*2), 2) - mean(redundant_traces_v1(:,fs/2+1:fs), 2);
dev_responses = mean(deviant_traces_v1(:,preSeconds*fs:fs*2), 2) - mean(deviant_traces_v1(:,fs/2+1:fs), 2);

means = [mean(control_responses) mean(red_responses) mean(dev_responses)];

y_vals = nan(3, length(red_responses));
y_vals(1, 1:length(control_responses)) = control_responses;
y_vals(2, 1:length(red_responses)) = red_responses;
y_vals(3, 1:length(dev_responses)) = dev_responses;

x_labels = categorical({'Control', 'Redundant', 'Deviant'}); 
x_labels = reordercats(x_labels, {'Control' 'Redundant' 'Deviant'});  

errorplus=nan(1, size(y_vals, 1));
for i = 1:size(y_vals, 1)
    errorplus(1, i) = std(y_vals(i,:), 'omitnan')/(sqrt(sum(~isnan(y_vals(i,:)))));
end
errorplus = reshape(errorplus, length(errorplus)/length(x_labels), length(x_labels));

% make fig 
barplot_fig = figure(); hold on
b = bar(x_labels, means);
errorbar(x_labels, means, errorplus, 'k.');
box off
exportgraphics(barplot_fig, 'crd_barplot.pdf','ContentType','vector')

% make fig 
barplot_fig = figure(); hold on
boxchart(y_vals','Notch','on');
box off
exportgraphics(barplot_fig, 'crd_boxplot.pdf','ContentType','vector')


save("barboxplot_data.mat", "x_labels", "means", "y_vals")



%%

iti = oddballStimuliData.interStimIntervalS;
stim_dur = oddballStimuliData.stimDurationS;
preSamples = ceil((iti+stim_dur)*10*3.25);
postSamples = floor((iti+stim_dur)*10*0.9);
%alignedControl = squeeze(mean(averageMovieToEvent(zData,mesoTimestamps,[controlDeviant1TimestampsOn; controlDeviant2TimestampsOn],preSamples,postSamples),1, 'omitnan'));
alignedDeviant = squeeze(mean(averageMovieToEvent(zData,mesoTimestamps,[stimB_dev_TimestampsOn; stimA_dev_TimestampsOn],preSamples,postSamples),1, 'omitnan'));
%alignedRedundant = squeeze(mean(averageMovieToEvent(zData,mesoTimestamps,[stimA_red_TimestampsOn],preSamples,postSamples),1, 'omitnan'));

dev_time = preSamples+1;
red_times = [preSamples-1*((iti+stim_dur)*10) preSamples-2*((iti+stim_dur)*10) preSamples-3*((iti+stim_dur)*10)];

%% align deviant response to v1

%align to deviantTimestamp
v1_deviant =  mean(alignedDeviant([1;2], :), 1, 'omitnan');


%plot
figure();
plot(v1_deviant);
xline(dev_time, 'r')
xline(red_times, 'g')
title('average v1 response')


%% lets look at all the facemap stuff 


% original full face pcs
figure(); hold on
for i = 1:20
    pc = face_proc.motSVD{1, 2}(:,i);
    if skewness(pc) < 0.5; pc = -pc; end
    face_Norm = normalize(pc(1:maxlength));
    
    subplot(5, 4, i)
    aligned_face = nan(length(stimB_dev_TimestampsOn), preSamples+postSamples, 'single');
    for ii = 1:length(stimB_dev_TimestampsOn)
        ind = findInSorted(pupil_time,stimB_dev_TimestampsOn(ii));
        if ind-preSamples> 0 && ind+postSamples-1 <= length(pupil_time)
            aligned_face(ii,:) = face_Norm(ind-preSamples:ind+postSamples-1);
        end
    end
    plot(squeeze(mean(aligned_face,1)), 'k')

end



% additional pc 
load(fullfile(session, "additional_face.mat"));
add_proc = proc;
figure(); hold on
for i = 1:20
    pc = add_proc.motSVD{1, 4}(:,i);
    if skewness(pc) < 0.5; pc = -pc; end
    face_Norm = normalize(pc(1:maxlength));
    
    subplot(5, 4, i)
    aligned_face = nan(length(stimB_dev_TimestampsOn), preSamples+postSamples, 'single');
    for ii = 1:length(stimB_dev_TimestampsOn)
        ind = findInSorted(pupil_time,stimB_dev_TimestampsOn(ii));
        if ind-preSamples> 0 && ind+postSamples-1 <= length(pupil_time)
            aligned_face(ii,:) = face_Norm(ind-preSamples:ind+postSamples-1);
        end
    end
    plot(squeeze(mean(aligned_face,1)), 'k')

end



%% let's make a r-r-r-d plot with all variables shown... 



figure(); hold on

subplot(3,1, 1)
plot(v1_deviant, 'k');
ylabel('dff');

subplot(3,1,2)
vis_trace = zeros((preSamples+postSamples), 1);
vis_trace(red_times(1):red_times(1)+stim_dur*10) = 1;
vis_trace(red_times(2):red_times(2)+stim_dur*10) = 1;
vis_trace(red_times(3):red_times(3)+stim_dur*10) = 1;
vis_trace(dev_time:dev_time+stim_dur*10) = 1;
plot(vis_trace, 'k')
ylabel('vis stim'); ylim([0 1.5])


subplot(3, 1, 3)
aligned_face = nan(length(stimB_dev_TimestampsOn), preSamples+postSamples, 'single');  
for i = 1:length(stimB_dev_TimestampsOn)
    ind = findInSorted(pupil_time,stimB_dev_TimestampsOn(i));
    if ind-preSamples> 0 && ind+postSamples-1 <= length(pupil_time)
        aligned_face(i,:) = face_Norm(ind-preSamples:ind+postSamples-1);
    end
end
plot(squeeze(mean(aligned_face,1)), 'k')
ylabel('FaceMap PC1')


%% make video 

d = alignedDeviant;
d = reshape(d, 256, 256, []);
for i=1:size(d, 3)
    temp = d(:,:,i);
    temp(boundaryMask) = nan;
    d(:,:,i) = temp;
end


filename = 'redundant_followedby_deviant.gif';

for i =1:size(d,3)
    h=figure('Position', [100 100 900 700]);
    % cortex wide vid
    toShow = d(:,:,i);
    imagesc(toShow)
    axis off
    caxis([-1 1]);
    % plot parcels
    title('redundant followed by deviant','FontSize',20)
    drawnow;
    frame = getframe(h); 
    im = frame2im(frame); 
    whiteMask = im == mode(im,[1 2]);
    tempIm = permute(im,[3 1 2]);
    tempIm(1,whiteMask(:,:,1)) = squeeze(im(1,1,1));
    tempIm(2,whiteMask(:,:,2)) = squeeze(im(1,1,2));
    tempIm(3,whiteMask(:,:,3)) = squeeze(im(1,1,3));
    whiteIm = permute(tempIm,[2 3 1]);
    [imind,cm] = rgb2ind(whiteIm,256); 
    if i == 1 
          imwrite(imind,cm,filename,'gif', 'Loopcount',inf,'DelayTime',.1); 
    else 
          imwrite(imind,cm,filename,'gif','WriteMode','append','DelayTime',.1); 
    end 
    close(h)
end









