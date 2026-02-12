
addpath(genpath('D:\meso_demo\analysis_scripts'))
addpath(genpath('D:\meso_demo\MesoProcessing-master'))


session = 'C:\Users\CardinLab\Desktop\Oddball_testing\59_platformtest\59_oddball_piezo_10132022\';
mkNewDir(fullfile(session, 'oddball'))
cd(fullfile(session, 'oddball'))

%% load neural signal 
load(fullfile(session, 'normed_lowface_blue_uv.mat'));
zData = (normed_sig-mean(normed_sig,2, 'omitnan'))./std(normed_sig,[],2, 'omitnan');

load('parcells_updated121519.mat');
boundaryMask = sum(parcells_new.indicators,3)>1;
v1_left=find(parcells_new.CombinedParcells==1);
v1_right=find(parcells_new.CombinedParcells==2);


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

if length(mergeOn) ~= 850
    warning('unexpected number of diode pulses')
    disp(['number of diode pulses: ' num2str(length(mergeOn))])
end

%% load facemap data
facemapfile = ls(fullfile(session, '\*proc.mat'));
load(fullfile(session, facemapfile), 'proc')

pc = proc.motSVD{1, 2}(:,1);
if skewness(pc) < 0.5; pc = -pc; end

%% define and trim data

maxlength = min(min(length(spike2_data.pupilFrameOnTimestamps), length(pc)), min(size(zData,2), length(spike2_data.blueOnTimestamps)));

pupil_time = spike2_data.pupilFrameOnTimestamps(1:maxlength);
face_Norm = normalize(pc(1:maxlength));
zData = zData(:, 1:maxlength);
mesoTimestamps = spike2_data.blueOnTimestamps(1:maxlength);


%% define face states

load(fullfile(session, 'state_timestamps.mat'))
sitOn = mesoTimestamps(1);
sitOff = mesoTimestamps(end);


% gets time points within sitting bouts when not grooming
[notGroomOn, notGroomOff] = timestampsNot(states.groomHighOn, states.groomHighOff, 'StartTime', sitOn, 'EndTime', sitOff);
[sitOn_final, sitOff_final] = timestampsAnd(sitOn, sitOff, notGroomOn, notGroomOff);

% gets time points within sitting bouts 
pupilTime_Idx=cell(1,length(sitOn_final));
for st=1:length(sitOn_final)
    pupilTime_Idx{st}=find(pupil_time > sitOn_final(st) & pupil_time < sitOff_final(st));
end
pupilTime=cell2mat(pupilTime_Idx');
face=face_Norm(pupilTime);

% define quantiles 
zThresh_bottom1 = quantile(face,0.75);
zThresh_bottom2 = quantile(face,0.5);
zThresh_bottom3 = quantile(face,0.25);




%% get oddball experiment parameters

load(fullfile(session, 'oddballStimuliData.mat'));
numResidualRedundants = oddballStimuliData.numPresentations-oddballStimuliData.numRedundantPrecedingDeviant*oddballStimuliData.numDeviants; % number ooddballStimuliData.spatialFreq redundant stimuli that aren't in the 3 stimuli beoddballStimuliData.spatialFreqore a deviant
manyControlsLength = oddballStimuliData.numPresentations;
epochLengths = [oddballStimuliData.numPresentations oddballStimuliData.preSequenceRedundantStimuli oddballStimuliData.numPresentations oddballStimuliData.preSequenceRedundantStimuli oddballStimuliData.numPresentations];
epochMasks = false(sum(epochLengths),length(epochLengths));

currStart = 0;
for ii = 1:length(epochLengths)
  epochMasks(currStart+1:currStart+epochLengths(ii),ii) = true;
  currStart = currStart + epochLengths(ii);
end

deviantStim1 = mode(oddballStimuliData.masterSequence(epochMasks(:,5)));
deviantStim2 = mode(oddballStimuliData.masterSequence(epochMasks(:,3)));


controlDeviant1TimestampsOn = mergeOn(epochMasks(:,1) & oddballStimuliData.masterSequence==deviantStim1);
controlDeviant1TimestampsOff = mergeOff(epochMasks(:,1) & oddballStimuliData.masterSequence==deviantStim1);

deviant1TimestampsOn = mergeOn(epochMasks(:,3) & oddballStimuliData.masterSequence==deviantStim1);
deviant1TimestampsOff = mergeOff(epochMasks(:,3) & oddballStimuliData.masterSequence==deviantStim1);
preStim1TimestampsOn = mergeOn(epochMasks(:,2));
preStim1TimestampsOff = mergeOff(epochMasks(:,2));

redundant1TimestampsOn = mergeOn(epochMasks(:,3) & oddballStimuliData.masterSequence==deviantStim2);
redundant1TimestampsOff = mergeOff(epochMasks(:,3) & oddballStimuliData.masterSequence==deviantStim2);

controlDeviant2TimestampsOn = mergeOn(epochMasks(:,1) & oddballStimuliData.masterSequence==deviantStim2);
controlDeviant2TimestampsOff = mergeOff(epochMasks(:,1) & oddballStimuliData.masterSequence==deviantStim2);
deviant2TimestampsOn = mergeOn(epochMasks(:,5) & oddballStimuliData.masterSequence==deviantStim2);
deviant2TimestampsOff = mergeOff(epochMasks(:,5) & oddballStimuliData.masterSequence==deviantStim2);
preStim2TimestampsOn = mergeOn(epochMasks(:,4));
preStim2TimestampsOff = mergeOff(epochMasks(:,4));

redundant2TimestampsOn = mergeOn(epochMasks(:,5) & oddballStimuliData.masterSequence==deviantStim1);
redundant2TimestampsOff = mergeOff(epochMasks(:,5) & oddballStimuliData.masterSequence==deviantStim1);


%% align deviant response to higher order visual areas

%first need to plot deviant response
avg_deviant = squeeze(mean(averageMovieToEvent(zData,mesoTimestamps,[deviant1TimestampsOn; deviant2TimestampsOn],0,1),1, 'omitnan'));
avg_deviant = squeeze(reshape(avg_deviant, [], 256, 256));
figure();
imagesc(avg_deviant);

% use these to draw box around rois.. then right click and select "create mask"
left_roi = roipoly;
right_roi = roipoly;

%align to deviantTimestamp
v2_left_deviant = mean(alignedDeviant(left_roi, :), 1, 'omitnan');
v2_right_deviant = mean(alignedDeviant(right_roi, :), 1, 'omitnan');
v2_deviant = [v2_left_deviant; v2_right_deviant];
v2_deviant = mean(v2_deviant, 1);

%plot
figure();
plot(v2_deviant);
title('average response of outlined rois')


%% split based on facemap 

% just to see what is going on... 
v1_left_zData = zData(v1_left, :);
v1_right_zData = zData(v1_right, :);
v1_data = mean([v1_left_zData;v1_right_zData], 'omitnan');

% figure(); hold on
% plot(mesoTimestamps, v1_data)
% plot([deviant1TimestampsOn; deviant2TimestampsOn], 6, '*g')

%%
redundant_Ts = [redundant1TimestampsOn; ];%redundant2TimestampsOn; 
redundant_amp_face = nan(length(redundant_Ts), 2);
for i = 1:length(redundant_Ts)

    ind = findInSorted(mesoTimestamps,redundant_Ts(i));
    for groombout = 1:length(states.groomHighOn)
        if mesoTimestamps(ind) > states.groomHighOn(groombout) && mesoTimestamps(ind) < states.groomHighOff(groombout)
            disp('true');
            continue;
        end
    end
    v1_amplitude_before = sum(v1_data(ind-9:ind));
    v1_amplitude_after = sum(v1_data(ind:ind+9));
    v1_amplitude = v1_amplitude_after - v1_amplitude_before;
    face_val = mean(face_Norm(ind-4:ind));
    redundant_amp_face(i, :) = [v1_amplitude, face_val];
    
end

redundant_amps_highface = redundant_amp_face(redundant_amp_face(:,2) > zThresh_bottom1, :);
redundant_amps_midhighface = redundant_amp_face(redundant_amp_face(:,2) < zThresh_bottom1 & redundant_amp_face(:,2) > zThresh_bottom2, :);
redundant_amps_midlowface = redundant_amp_face(redundant_amp_face(:,2) < zThresh_bottom2 & redundant_amp_face(:,2) > zThresh_bottom3, :);
redundant_amps_lowface = redundant_amp_face(redundant_amp_face(:,2) < zThresh_bottom3, :);

redundant_amps.redundant_amps_highface=redundant_amps_highface;
redundant_amps.redundant_amps_midhighface=redundant_amps_midhighface;
redundant_amps.redundant_amps_midlowface=redundant_amps_midlowface;
redundant_amps.redundant_amps_lowface=redundant_amps_lowface;
save('redundant_amps.mat', "redundant_amps");


% regression figs 
redunant_reg = figure();
tiledlayout(2,2, 'TileSpacing', 'tight', 'Padding', 'compact');

all_vals = [redundant_amps_lowface(:,1); redundant_amps_midlowface(:,1); redundant_amps_midhighface(:,1); redundant_amps_highface(:,1)];
minVal = min(all_vals);
maxVal = max(all_vals);

nexttile(1);
scatterPlot(redundant_amps_lowface(:,1), minVal, maxVal)
title('AUC 1s following redundant stim, low face')

nexttile(2);
scatterPlot(redundant_amps_midlowface(:,1), minVal, maxVal)
title('AUC 1s following redundant stim, low-mid face')

nexttile(3);
scatterPlot(redundant_amps_midhighface(:,1), minVal, maxVal)
title('AUC 1s following redundant stim, high-mid face')

nexttile(4);
scatterPlot(redundant_amps_highface(:,1), minVal, maxVal)
title('AUC 1s following redundant stim, high face')

saveas(redunant_reg, 'redundant_regs_fig.fig')



%bar plot 
x=categorical({'Low face', 'Low-mid face', 'High-mid face', 'High face'});
x=reordercats(x, {'Low face', 'Low-mid face', 'High-mid face', 'High face'});
y_means = [mean(redundant_amps_midlowface(:,1)); mean(redundant_amps_midlowface(:,1)); mean(redundant_amps_midhighface(:,1)); mean(redundant_amps_highface(:,1)); ];
barplot = figure; hold on; box off
bar(x,y_means); 
ylabel('AUC 1s following vis stim') %,'FontSize',18
title('redundant')

saveas(barplot, 'barplot_AUC_redundant.png')




%% now looking at deviant

deviant_Ts = [deviant1TimestampsOn; ];%deviant2TimestampsOn
deviant_amp_face = nan(length(deviant_Ts), 2);
for i = 1:length(deviant_Ts)

    ind = findInSorted(mesoTimestamps,deviant_Ts(i));
    for groombout = 1:length(states.groomHighOn)
        if  mesoTimestamps(ind) > states.groomHighOn(groombout) &&  mesoTimestamps(ind) < states.groomHighOff(groombout)
            disp('true');
            continue;
        end
    end
    v1_amplitude_before = sum(v1_data(ind-9:ind));
    v1_amplitude_after = sum(v1_data(ind:ind+9));
    v1_amplitude = v1_amplitude_after - v1_amplitude_before;
    face_val = mean(face_Norm(ind-4:ind));
    deviant_amp_face(i, :) = [v1_amplitude, face_val];
    
end

deviant_amps_highface = deviant_amp_face(deviant_amp_face(:,2) > zThresh_bottom1, :);
deviant_amps_midhighface = deviant_amp_face(deviant_amp_face(:,2) < zThresh_bottom1 & deviant_amp_face(:,2) > zThresh_bottom2, :);
deviant_amps_midlowface = deviant_amp_face(deviant_amp_face(:,2) < zThresh_bottom2 & deviant_amp_face(:,2) > zThresh_bottom3, :);
deviant_amps_lowface = deviant_amp_face(deviant_amp_face(:,2) < zThresh_bottom3, :);

deviant_amps.redundant_amps_highface=deviant_amps_highface;
deviant_amps.redundant_amps_midhighface=deviant_amps_midhighface;
deviant_amps.redundant_amps_midlowface=deviant_amps_midlowface;
deviant_amps.redundant_amps_lowface=deviant_amps_lowface;
save('deviant_amps.mat', "deviant_amps");


% regression figs 
deviant_reg = figure();
tiledlayout(2,2, 'TileSpacing', 'tight', 'Padding', 'compact');

all_vals = [deviant_amps_lowface(:,1); deviant_amps_midlowface(:,1); deviant_amps_midhighface(:,1); deviant_amps_highface(:,1)];
minVal = min(all_vals);
maxVal = max(all_vals);

nexttile(1);
scatterPlot(deviant_amps_lowface(:,1), minVal, maxVal)
title('AUC 1s following deviant stim, low face')

nexttile(2);
scatterPlot(deviant_amps_midlowface(:,1), minVal, maxVal)
title('AUC 1s following deviant stim, low-mid face')

nexttile(3);
scatterPlot(deviant_amps_midhighface(:,1), minVal, maxVal)
title('AUC 1s following deviant stim, high-mid face')

nexttile(4);
scatterPlot(deviant_amps_highface(:,1), minVal, maxVal)
title('AUC 1s following deviant stim, high face')

saveas(deviant_reg, 'deviant_regs_fig.fig')


%bar plot 
x=categorical({'Low face', 'Low-mid face', 'High-mid face', 'High face'});
x=reordercats(x, {'Low face', 'Low-mid face', 'High-mid face', 'High face'});
y_means = [mean(deviant_amps_midlowface(:,1)); mean(deviant_amps_midlowface(:,1)); mean(deviant_amps_midhighface(:,1)); mean(deviant_amps_highface(:,1)); ];
barplot = figure; hold on; box off
bar(x,y_means); 
ylabel('AUC 1s following vis stim') %,'FontSize',18
title('deviant')


saveas(barplot, 'barplot_AUC_deviant_firstblockonly.png')


%% make bar plot like in hamm et al

control_Ts = [controlDeviant1TimestampsOn; controlDeviant2TimestampsOn];
control_amp_face = nan(length(control_Ts), 2);
for i = 1:length(control_Ts)

    ind = findInSorted(mesoTimestamps,control_Ts(i));
    for groombout = 1:length(states.groomHighOn)
        if  mesoTimestamps(ind) > states.groomHighOn(groombout) &&  mesoTimestamps(ind) < states.groomHighOff(groombout)
            disp('true');
            continue;
        end
    end
    v1_amplitude_before = sum(v1_data(ind-9:ind));
    v1_amplitude_after = sum(v1_data(ind:ind+9));
    v1_amplitude = v1_amplitude_after - v1_amplitude_before;
    face_val = mean(face_Norm(ind-4:ind));
    control_amp_face(i, :) = [v1_amplitude, face_val];
    
end

%bar plot 
x=categorical({'Control', 'Redundant', 'Deviant'});
x=reordercats(x, {'Control', 'Redundant', 'Deviant'});
y_means = [mean(control_amp_face(:,1)); mean(redundant_amp_face(:,1)); mean(deviant_amp_face(:,1))];
barplot = figure; hold on; box off
bar(x,y_means); 
ylabel('AUC 1s following vis stim') %,'FontSize',18
saveas(barplot, 'CvsRvsD.png')


%% lets look at piezo


time = spike2_data.analog_signal_time_vect;
hands = proc.motSVD{1, 3}(:,1);
if skewness(hands) < 0.5;  hands = -hands; end


cutoff_freq = 10;
sampling_freq = round(1/median(diff(time)));
[b,a] = butter(2,cutoff_freq/(sampling_freq/2), "low");
piezo = filtfilt(b,a,double(spike2_data.piezo));
%piezo = bandpass(spike2_data.piezo, [0.5 1.5], 5000);
piezo = abs(hilbert(piezo));
piezo = smooth(piezo, 500);
piezo = normalize(piezo);
figure(); plot(piezo)


% figure();
% subplot(2,1,1);
% plot(spike2_data.pupilFrameOnTimestamps(1:length(hands)), hands); hold on;
% xline(states.groomHighOn, 'g');
% xline(states.groomHighOff, 'r');
% subplot(2,1,2)
% plot(piezo); hold on
% plot(redundant_Ts*5000, 0.2, 'g*')
% plot(deviant_Ts*5000, 0.2, 'r*')


piezo_trace = figure();
ax1 = subplot(2,1,1);
plot(spike2_data.piezo)
ax2 = subplot(2,1,2);
plot(piezo)
linkaxes([ax1 ax2], 'x')
saveas(piezo_trace, 'piezo_trace.fig')


deviant_Ts = [deviant1TimestampsOn; deviant2TimestampsOn]; %
deviant_vidgets = nan(length(deviant_Ts), 1);
for i = 1:length(deviant_Ts)
    ind = findInSorted(time,deviant_Ts(i));
    for groombout = 1:length(states.groomHighOn)
        if  time(ind) > states.groomHighOn(groombout) &&  time(ind) < states.groomHighOff(groombout)
            disp('true');
            continue;
        end
    end
    piezo_val = sum(piezo(ind:ind+10000));
    deviant_vidgets(i) = piezo_val;
end

redundant_Ts = [redundant2TimestampsOn; redundant1TimestampsOn]; %
redundant_videgts = nan(length(redundant_Ts), 1);
for i = 1:length(redundant_Ts)
    ind = findInSorted(time,redundant_Ts(i));
    for groombout = 1:length(states.groomHighOn)
        if  time(ind) > states.groomHighOn(groombout) &&  time(ind) < states.groomHighOff(groombout)
            disp('true');
            continue;
        end
    end
    piezo_val = sum(piezo(ind:ind+10000));
    redundant_videgts(i) = piezo_val;
end

x=categorical({'Redundant', 'Deviant'});
x=reordercats(x, {'Redundant', 'Deviant'});
y_means = [mean(redundant_videgts); mean(deviant_vidgets)];
barplot = figure; hold on; box off
bar(x,y_means); 
title('integral of 2s following stim onset')
saveas(barplot, '2s_integral_following_onset.png')


redundant_vidgets = nan(length(redundant_Ts), 15000);
for i = 1:length(redundant_Ts)
    ind = findInSorted(time,redundant_Ts(i));
    for groombout = 1:length(states.groomHighOn)
        if  time(ind) > states.groomHighOn(groombout) &&  time(ind) < states.groomHighOff(groombout)
            disp('true');
            continue;
        end
    end
    piezo_trace = piezo(ind-4999:ind+10000);
    redundant_vidgets(i, :) = piezo_trace;
end

mean_redundant_vidget=mean(redundant_vidgets, 1, 'omitnan');
redundant_vidget_fig = figure();
plot(mean_redundant_vidget)
title('mean redundant vidget')
saveas(redundant_vidget_fig, 'redundant_vidget_fig.png')


deviant_vidgets = nan(length(deviant_Ts), 15000);
for i = 1:length(deviant_Ts)
    ind = findInSorted(time,deviant_Ts(i));
    piezo_trace = piezo(ind-4999:ind+10000);
    deviant_vidgets(i, :) = piezo_trace;
end

mean_deviant_vidget=mean(deviant_vidgets, 1, 'omitnan');
deviant_vidget_fig = figure();
plot(mean_deviant_vidget)
title('mean deviant vidget')
saveas(deviant_vidget_fig, 'deviant_vidget_fig.png')



%%

figure(); hold on
for i = 1:35
    subplot(7, 5, i)
    plot(redundant_vidgets(i, :))
    ylim([-0.5 15])
end




figure(); hold on
for i = 1:35
    subplot(7, 5, i)
    plot(deviant_vidgets(i, :))
    ylim([-0.5 5])
end


figure(); hold on
for i = 36:70
    subplot(7, 5, i-35)
    plot(deviant_vidgets(i, :))
    ylim([-0.5 5])
end








%%

iti = oddballStimuliData.interStimIntervalS;
stim_dur = oddballStimuliData.stimDurationS;
preSamples = ceil((iti+stim_dur)*10*.5);
postSamples = floor((iti+stim_dur)*10*.5);
%alignedControl = squeeze(mean(averageMovieToEvent(zData,mesoTimestamps,[controlDeviant1TimestampsOn; controlDeviant2TimestampsOn],preSamples,postSamples),1, 'omitnan'));
alignedDeviant = squeeze(mean(averageMovieToEvent(zData,mesoTimestamps,[deviant1TimestampsOn],preSamples,postSamples),1, 'omitnan'));
alignedRedundant = squeeze(mean(averageMovieToEvent(zData,mesoTimestamps,[redundant2TimestampsOn],preSamples,postSamples),1, 'omitnan'));



%% align deviant response to v1

%align to deviantTimestamp
v1_left_deviant = mean(alignedDeviant(v1_left, :), 1, 'omitnan');
v1_right_deviant = mean(alignedDeviant(v1_right, :), 1, 'omitnan');
v1_deviant = [v1_left_deviant; v1_right_deviant];
v1_deviant = mean(v1_deviant, 1);

%plot
figure();
plot(v1_deviant);
title('average v1 response')

figure();
b=reshape(alignedRedundant, 256, 256, []);
imagesc(b(:,:,33))


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



%%
d = alignedDeviant-alignedRedundant;
d = reshape(d, 256, 256, []);
for i=1:size(d, 3)
    temp = d(:,:,i);
    temp(boundaryMask) = nan;
    d(:,:,i) = temp;
end


filename = 'deviantSubRedundant.gif';

for i =1:size(d,3)
    h=figure('Position', [100 100 900 700]);
    % cortex wide vid
    toShow = d(:,:,i);
    imagesc(toShow)
    axis off
    caxis([-1 1]);
    % plot parcels
    title('Deviant-Redundant','FontSize',20)
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

%%

d = alignedDeviant-alignedControl;
d = reshape(d, 256, 256, []);
for i=1:30
    temp = d(:,:,i);
    temp(boundaryMask) = nan;
    d(:,:,i) = temp;
end




filename = 'deviantSubControl_try2.gif';

for i =1:size(d,3)
    h=figure('Position', [100 100 1500 700]);
    % cortex wide vid
    toShow = d(:,:,i);
    imagesc(toShow)
    axis off
    caxis([-1 1]);
    % plot parcels
    title('Deviant-Control','FontSize',20)
    drawnow;
    frame = getframe(h); 
    im = frame2im(frame); 
    [imind,cm] = rgb2ind(im,256); 
    if i == 1 
          imwrite(imind,cm,filename,'gif', 'Loopcount',inf,'DelayTime',.1); 
    else 
          imwrite(imind,cm,filename,'gif','WriteMode','append','DelayTime',.1); 
    end 
    close(h)
end










