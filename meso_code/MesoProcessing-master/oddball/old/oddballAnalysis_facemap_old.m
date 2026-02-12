
addpath(genpath('D:\meso_demo\analysis_scripts'))
addpath(genpath('D:\meso_demo\MesoProcessing-master'))
load('parcells_updated121519.mat');
boundaryMask = sum(parcells_new.indicators,3)>1;


session = 'D:\meso_demo\Oddball_testing\59_platformtest\59_oddball_09292022';
%session = 'D:\meso_demo\CDKL5_oddball\ro110-263\oddball_263_09142022\'; %cdkl5(+/-) ---> best looking one by far 
%session = 'D:\meso_demo\CDKL5_oddball\ro110-264\264_oddball_09142022\'; %cdkl5(-/-)


%% load neural signal 
load(fullfile(session, 'normed_lowface_blue_uv.mat'));
zData = (normed_sig-mean(normed_sig,2, 'omitnan'))./std(normed_sig,[],2, 'omitnan');

%% load and process spike2 
load(fullfile(session, 'final_timestamps.mat'));
[mergeOn,mergeOff] = mergeTimestamps(spike2_data.diodeOnTimestamps,spike2_data.diodeOffTimestamps,1);
if length(mergeOn) ~= 850
    warning('unexpected number of diode pulses')
    disp(['number of diode pulses: ' num2str(length(mergeOn))])
end

%% load facemap data
facemapfile = ls(fullfile(session, '\*proc.mat'));
load(fullfile(session, facemapfile), 'proc')

pc = proc.motSVD{1, 2}(:,1);
if skewness(pc) < 0.5; pc = -pc; end

% define and trim
imaging_time = spike2_data.blueOnTimestamps;
maxlength = min(length(spike2_data.pupilFrameOnTimestamps), length(pc));
pupil_time = spike2_data.pupilFrameOnTimestamps(1:maxlength);
face_Norm = normalize(pc);

%% define face states

load(fullfile(session, 'state_timestamps.mat'))
sitOn = imaging_time(1);
sitOff = imaging_time(end);
minStateDur = 2; 

% gets time points within sitting bouts when not grooming
[notGroomOn, notGroomOff] = timestampsNot(states.groomHighOn, states.groomHighOff, 'StartTime', sitOn, 'EndTime', sitOff);
[sitOn_final, sitOff_final] = timestampsAnd(sitOn, sitOff, notGroomOn, notGroomOff);

% gets time points within sitting bouts 
pupilTime_Idx=cell(1,length(sitOn_final));
for st=1:length(sitOn_final)
    pupilTime_Idx{st}=find(pupil_time > sitOn_final(st) & pupil_time < sitOff_final(st));
end
pupilTime_quiescence=cell2mat(pupilTime_Idx');
face_quiescence=face_Norm(pupilTime_quiescence);

% define quantiles 
zThresh_top = quantile(face_quiescence,0.9);
zThresh_bottom1 = quantile(face_quiescence,0.7);
zThresh_bottom2 = quantile(face_quiescence,0.5);
zThresh_bottom3 = quantile(face_quiescence,0.3);
zThresh_bottom4 = quantile(face_quiescence,0.1);


%% maybe instead of this... just assign values between thresholds and if vis stim happens within that then assign accordingly


% Top 25% 
[face25on,face25off] = cjb_changepoints_bw(face_Norm, pupil_time, zThresh_top, zThresh_bottom1, 1); 
[face25on_inter,face25off_inter] = timestampsAnd(face25on,face25off,sitOn,sitOff);
face25onFinal=squeeze(face25on_inter(face25off_inter - face25on_inter > minStateDur));
face25offFinal=squeeze(face25off_inter(face25off_inter - face25on_inter > minStateDur));

% 25-50%
[face50on,face50off] = cjb_changepoints_bw(face_Norm, pupil_time, zThresh_bottom1, zThresh_bottom2, 1); 
[face50on_inter,face50off_inter] = timestampsAnd(face50on,face50off,sitOn,sitOff);
face50onFinal=squeeze(face50on_inter(face50off_inter - face50on_inter > minStateDur));
face50offFinal=squeeze(face50off_inter(face50off_inter - face50on_inter > minStateDur));

% 50-75%
[face75on,face75off] = cjb_changepoints_bw(face_Norm, pupil_time, zThresh_bottom2, zThresh_bottom3, 1); 
[face75on_inter,face75off_inter] = timestampsAnd(face75on,face75off,sitOn,sitOff);
face75onFinal=squeeze(face75on_inter(face75off_inter - face75on_inter > minStateDur));
face75offFinal=squeeze(face75off_inter(face75off_inter - face75on_inter > minStateDur));

% bottom 25% 
[face100on,face100off] = cjb_changepoints_bw(face_Norm, pupil_time, zThresh_bottom3, zThresh_bottom4, 1); 
[face100on_inter,face100off_inter] = timestampsAnd(face100on,face100off,sitOn,sitOff);
face100onFinal=squeeze(face100on_inter(face100off_inter - face100on_inter > minStateDur));
face100offFinal=squeeze(face100off_inter(face100off_inter - face100on_inter > minStateDur));


figure();
[b,a] = butter(2,0.1/(10/2));%c cutoff_freq/(samplingfreq/2)
pc_filt = filtfilt(b,a,double(face_Norm));

ax(1) = subplot(4, 1, 1);
plot(spike2_data.pupilFrameOnTimestamps(1:maxlength), face_Norm(1:maxlength)); hold on
xline(face25onFinal, 'g')
xline(face25offFinal, 'r')

ax(2) = subplot(4, 1, 2);
plot(spike2_data.pupilFrameOnTimestamps(1:maxlength), face_Norm(1:maxlength)); hold on
xline(face50onFinal, 'g')
xline(face50offFinal, 'r')

ax(3) = subplot(4, 1, 3);
plot(spike2_data.pupilFrameOnTimestamps(1:maxlength), face_Norm(1:maxlength)); hold on
xline(face75onFinal, 'g')
xline(face75offFinal, 'r')

ax(4) = subplot(4, 1, 4);
plot(spike2_data.pupilFrameOnTimestamps(1:maxlength), face_Norm(1:maxlength)); hold on
xline(face100onFinal, 'g')
xline(face100offFinal, 'r')

linkaxes(ax,'x');

%[length(face25offFinal), length(face50offFinal), length(face75offFinal), length(face100offFinal)]



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




%%
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



%%
maxLength = min(length(spike2_data.pupilFrameOffTimestamps), length(spike2_data.blueOnTimestamps));
mesoTimestamps = spike2_data.blueOnTimestamps(1:maxLength);
iti = oddballStimuliData.interStimIntervalS;
stim_dur = oddballStimuliData.stimDurationS;
preSamples = (iti+stim_dur)*10*3;
postSamples = (iti+stim_dur)*10*1;
%alignedControl = squeeze(mean(averageMovieToEvent(zData,mesoTimestamps,[controlDeviant1TimestampsOn; controlDeviant2TimestampsOn],preSamples,postSamples),1, 'omitnan'));
alignedDeviant = squeeze(mean(averageMovieToEvent(zData,mesoTimestamps,[deviant1TimestampsOn; deviant2TimestampsOn],preSamples,postSamples),1, 'omitnan'));
%alignedRedundant = squeeze(mean(averageMovieToEvent(zData,mesoTimestamps,[redundant1TimestampsOn; redundant2TimestampsOn],preSamples,postSamples),1, 'omitnan'));



%%

% v1 and v2 
v1_left=find(parcells_new.CombinedParcells==1);
v1_right=find(parcells_new.CombinedParcells==2);

%align to deviantTimestamp
v1_left_deviant = mean(alignedDeviant(v1_left, :), 1, 'omitnan');
v1_right_deviant = mean(alignedDeviant(v1_right, :), 1, 'omitnan');
v1_deviant = [v1_left_deviant; v1_right_deviant];
v1_deviant = mean(v1_deviant, 1);

%plot
figure();
plot(v1_deviant);
title('59 09292022 platform')

%split based on facemap 






%%
d = alignedDeviant;
d = reshape(d, 256, 256, []);
for i=1:30
    temp = d(:,:,i);
    temp(boundaryMask) = nan;
    d(:,:,i) = temp;
end


filename = 'deviant_try2.gif';

for i =1:size(d,3)
    h=figure('Position', [100 100 900 700]);
    % cortex wide vid
    toShow = d(:,:,i);
    imagesc(toShow)
    axis off
    caxis([-1 1]);
    % plot parcels
    title('Deviant','FontSize',20)
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
for i=1:30
    temp = d(:,:,i);
    temp(boundaryMask) = nan;
    d(:,:,i) = temp;
end


filename = 'deviantSubRedundant_try2.gif';

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










