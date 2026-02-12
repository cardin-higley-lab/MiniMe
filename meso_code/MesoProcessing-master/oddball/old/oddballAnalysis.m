
addpath(genpath('D:\meso_demo\analysis_scripts'))
addpath(genpath('D:\meso_demo\MesoProcessing-master'))
load('parcells_updated121519.mat');

session = 'D:\meso_demo\oddball\62_platformtest\oddball_platform_62_09282022\';

%% load neural signal 
load(fullfile(session, 'normed_lowface_blue_uv.mat'));
zData = (normed_sig-nanmean(normed_sig,3))./nanstd(normed_sig,[],3);

%% load and process spike2 and facemap data
load(fullfile(session, 'final_timestamps.mat'));
[mergeOn,mergeOff] = mergeTimestamps(spike2_data.diodeOnTimestamps,spike2_data.diodeOffTimestamps,1);
if length(mergeOn) ~= 850
    warning('unexpected number of diode pulses')
    disp(['number of diode pulses: ' num2str(length(mergeOn))])
end

facemapfile = ls(fullfile(session, '\*proc.mat'));
load(fullfile(session, facemapfile), 'proc')

pc = proc.motSVD{1, 2}(:,1);
if skewness(pc) < 0.5; pc = -pc; end

% define and trim
imaging_time = spike2_data.blueOnTimestamps;
maxlength = min(length(spike2_data.pupilFrameOnTimestamps), length(pc));
pupil_time = spike2_data.pupilFrameOnTimestamps(1:maxlength);
face_Norm = normalize(pc);

cutoff_freq = 0.1;
sampling_freq = round(1/median(diff(spike2_data.pupilFrameOnTimestamps)));
[b,a] = butter(2,cutoff_freq/(sampling_freq/2));
pc_filt = filtfilt(b,a,double(pc));

%% define face states

minStateDur = 2; 
load(fullfile(session, 'state_timestamps.mat'))
groomHighOn = states.groomHighOn;
groomHighOff = states.groomHighOff;
sitOn = imaging_time(1);
sitOff = imaging_time(end);


%% FACE MAP

% gets time points within sitting bouts when not grooming

[notGroomOn, notGroomOff] = timestampsNot(groomHighOn, groomHighOff, 'StartTime', sitOn, 'EndTime', sitOff);
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
zThresh_bottom1=quantile(face_quiescence,0.7);
zThresh_bottom2=quantile(face_quiescence,0.5);
zThresh_bottom3=quantile(face_quiescence,0.3);
zThresh_bottom4 = quantile(face_quiescence,0.1);


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
ax(1) = subplot(4, 1, 1);
plot(spike2_data.pupilFrameOnTimestamps(1:maxlength), pc_filt(1:maxlength)); hold on
xline(face25onFinal, 'g')
xline(face25offFinal, 'r')

ax(2) = subplot(4, 1, 2);
plot(spike2_data.pupilFrameOnTimestamps(1:maxlength), pc_filt(1:maxlength)); hold on
xline(face50onFinal, 'g')
xline(face50offFinal, 'r')

ax(3) = subplot(4, 1, 3);
plot(spike2_data.pupilFrameOnTimestamps(1:maxlength), pc_filt(1:maxlength)); hold on
xline(face75onFinal, 'g')
xline(face75offFinal, 'r')

ax(4) = subplot(4, 1, 4);
plot(spike2_data.pupilFrameOnTimestamps(1:maxlength), pc_filt(1:maxlength)); hold on
xline(face100onFinal, 'g')
xline(face100offFinal, 'r')

linkaxes(ax,'x');




[length(face25offFinal), length(face50offFinal), length(face75offFinal), length(face100offFinal)]



































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
mesoTimestamps = spike2_data.blueOnTimestamps(1:size(zData,3));
alignedControl = squeeze(nanmean(averageMovieToEvent(zData,mesoTimestamps,[controlDeviant1TimestampsOn; controlDeviant2TimestampsOn],10,20),1));
alignedDeviant = squeeze(nanmean(averageMovieToEvent(zData,mesoTimestamps,[deviant1TimestampsOn; deviant2TimestampsOn],10,20),1));
alignedRedundant = squeeze(nanmean(averageMovieToEvent(zData,mesoTimestamps,[redundant1TimestampsOn; redundant2TimestampsOn],10,20),1));




d = alignedDeviant-alignedRedundant;
boundaryMask = sum(parcells_new.indicators,3)>1;
for i=1:30
    temp = d(:,:,i);
    temp(boundaryMask) = nan;
    d(:,:,i) = temp;
end


filename = 'deviantSubRedundant.gif';

for i =1:size(d,3)
    h=figure('Position', [100 100 1500 700]);
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
boundaryMask = sum(parcells_new.indicators,3)>1;
for i=1:30
    temp = d(:,:,i);
    temp(boundaryMask) = nan;
    d(:,:,i) = temp;
end




filename = 'deviantSubControl.gif';

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
%%
tiledlayout(2,1);
ax(1) = nexttile;
plot(mesoTimestamps,trace);
ax(2) = nexttile;
plot((0:length(spike2_data.wheelSpeed)-1)/5000,spike2_data.wheelSpeed);
linkaxes(ax,'x')


%% prepping for symbolic regression

trace = squeeze(mean(zData(58:72,64:74,:),[1 2]));

stateVect = spike2_data.wheelSpeed>0.05;
stateVect(1) = 0;
wheelOnInds = find(diff(stateVect)==1)+1;
wheelOffInds = find(diff(stateVect)==-1)+1;
wheelOnInds = wheelOnInds(1:length(wheelOffInds));
speeds = zeros(size(wheelOnInds));
for i =1:length(speeds)
    speeds(i) = mean(spike2_data.wheelSpeed(wheelOnInds(i):wheelOffInds(i)));
end

toKeep = speeds>10;
toKeep2 = (wheelOffInds-wheelOnInds)>(3*5000);

t = (0:length(spike2_data.wheelSpeed)-1)/5000;

[wheelOn,wheelOff] = mergeTimestamps(t(wheelOnInds((toKeep&toKeep2))),t(wheelOffInds((toKeep&toKeep2))),1);
loco = zeros(size(trace));
for i =1:length(wheelOn)
    loco(mesoTimestamps>=wheelOn(i) & mesoTimestamps<wheelOff(i)) = 1;
end
loco2 = spike2_data.wheelSpeed(round(mesoTimestamps*5000));
%loco = pSpike2Data.wheelSpeed(round(mesoTimestamps*5000));

lowAbundanceTimestampsOn = [mergeOn(epochMasks(:,1)); deviant1TimestampsOn; deviant2TimestampsOn];
lowAbundanceTimestampsOff = [mergeOff(epochMasks(:,1)); deviant1TimestampsOff; deviant2TimestampsOff];
lowAbundanceStim = zeros(size(trace));
for i =1:length(lowAbundanceTimestampsOn)
    lowAbundanceStim(mesoTimestamps>=lowAbundanceTimestampsOn(i) & mesoTimestamps<lowAbundanceTimestampsOff(i)) = 1;
end

highAbundanceTimestampsOn = [preStim1TimestampsOn; preStim2TimestampsOn; redundant1TimestampsOn; redundant2TimestampsOn];
highAbundanceTimestampsOff = [preStim1TimestampsOff; preStim2TimestampsOff; redundant1TimestampsOff; redundant2TimestampsOff];
highAbundanceStim = zeros(size(trace));
for i =1:length(highAbundanceTimestampsOn)
    highAbundanceStim(mesoTimestamps>=highAbundanceTimestampsOn(i) & mesoTimestamps<highAbundanceTimestampsOff(i)) = 1;
end


X = zeros(length(trace),4);
X(:,1) = loco2;
X(:,2) = lowAbundanceStim;
X(:,3) = highAbundanceStim;
X(:,4) = trace;
X = X(1001:end,:);
save('X','X')


plot(mesoTimestamps(1001:end),pred_dot(:,4))
scatterTimestamps(highAbundanceTimestampsOn,highAbundanceTimestampsOff,5)
scatterTimestamps(lowAbundanceTimestampsOn,lowAbundanceTimestampsOff,6)

figure
preSamples = 10;
postSamples = 10;

alignedLowAbundance = squeeze(averageParcelsToEvent(pred_dot(:,4)',mesoTimestamps(1001:end),lowAbundanceTimestampsOn,preSamples,postSamples));
alignedHighAbundance = squeeze(averageParcelsToEvent(pred_dot(:,4)',mesoTimestamps(1001:end),highAbundanceTimestampsOn,preSamples,postSamples));

confidenceBand((-9:10)/10,alignedLowAbundance,0.95,1,'lineprops','g')
confidenceBand((-9:10)/10,alignedHighAbundance,0.95,1,'lineprops','r')
xline(0,'--k')
xline(.5,'--k')

figure

preSamples = 10;
postSamples = 20;

alignedLowAbundance = squeeze(averageParcelsToEvent(trace',mesoTimestamps,lowAbundanceTimestampsOn,preSamples,postSamples));
alignedHighAbundance = squeeze(averageParcelsToEvent(trace',mesoTimestamps,highAbundanceTimestampsOn,preSamples,postSamples));

confidenceBand((-10:19)/10,alignedLowAbundance,0.95,1,'lineprops','g')
confidenceBand((-10:19)/10,alignedHighAbundance,0.95,1,'lineprops','r')
xline(0,'--k')
xline(.5,'--k')

%%


dZ/dt = -0.283 + 3.7030*L + 2.8381*H + -0.809*Z  + 0.998*L*Z + 0.810*H*Z;

% no stim
dZ/dt = -0.809*Z -0.283

% Low Abundance Stim (control, oddball)
dZ/dt = -0.809*Z - 0.283  + 3.7030 + 0.998*Z = 0.1890*Z + 3.4200

% High Abundance Stim (redundant)
dZ/dt = -0.809*Z -0.283 + 2.8381 + 0.810*Z = 0.001*Z + 2.5551
