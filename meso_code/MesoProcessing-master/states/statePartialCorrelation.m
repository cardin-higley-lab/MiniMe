%% big iteration through data

% iterate through data in WorkingData
% get states with stateTimestamps
% parcel partialcorrelations by state
load('parcel_Obj.mat');




masterFolder = 'C:\Users\CardinLab\Desktop\PreprocessedData\Ketamine';
subFolders = dir([masterFolder '\*\*']);

added = 0; % iterator for added datasets (wont be same as folder, since ignore some datasets)
for folder = 1:length(subFolders)
    % ignore non data folders (.,..)
    if length(subFolders(folder).name)<3
        continue
    end
    disp(['Processing: ' subFolders(folder).name]);
    % get path to load data
    dataFolderPath = [subFolders(folder).folder '\' subFolders(folder).name];
    % load parcels, timestamps, and facemap data
    load(fullfile(dataFolderPath,'final_dFoF_parcels.mat'));
    facemapFileName=ls(fullfile(dataFolderPath, '\*proc.mat'));
    facemap = load(fullfile(dataFolderPath,facemapFileName));
    final_timestampsObj = load(fullfile(dataFolderPath,'final_timestamps.mat'));
    % older processing pipeline saved spike2_data as channel_data, so extract
    % accordingly
    if isfield(final_timestampsObj,'channel_data')
        spike2_data = final_timestampsObj.channel_data;
    else
        spike2_data = final_timestampsObj.spike2_data;
    end
    % facemap output from python and matlab is different. matlab output saves a
    if isfield(facemap,'proc')
        temp = facemap.proc.motSVD{end};
        pc = temp(:,1);
    else
        pc = facemap.motSVD_1(:,1);
    end
    % correct for sign of pc (depends on lighting conditions)
    if skewness(pc) < 0.5
        pc = -pc;
    end
 
    % Ensure no frame drops
    if abs(length(spike2_data.blueOnTimestamps)-length(dFoF_parcells.blue))>1 % give 1 as slack
        continue;
    elseif abs(length(spike2_data.pupilFrameOnTimestamps)-length(pc))>29 % drops last 27, give extra as slack
        continue;
    end
    
    % get timestamps for states
    states = stateTimestamps(spike2_data,pc);
    %load('parcel_Obj.mat');
    [numParcels,numSamples] = size(dFoF_parcells.blue);
    % some parcel traces will be nans, like colliculi, so find them
    nanlessParcelInds = ~isnan(dFoF_parcells.blue(:,1));
    %globalMean = zeros(numSamples,1);
    % get size of parcel so we can compute mean signal without loading
    % pixelwise data (takes forever)
    nanlessParcels = dFoF_parcells.blue(nanlessParcelInds,:);
    weightings = squeeze(sum(parcel_obj.indicators(:,:,nanlessParcelInds),[1 2])/sum(parcel_obj.indicators(:,:,nanlessParcelInds),'all'));
    meanTrace = sum(nanlessParcels.*weightings,1);
    % get non nan parcels to make things easier
    
    % take midpoint of exposure time as timepoinst
    mesoTimestamps = (spike2_data.blueOnTimestamps+spike2_data.blueOffTimestamps)/2; % really should use mesoFrameTimestamps, but seems to be having issues
    % dyanmically will grab data from states. useful if we choose to
    % increase the number of states we use (should require no modification
    % to this function)
    stateFields = fields(states);
    % have on and off timestamps, so /2 to get number of states
    numStates = length(stateFields)/2;
    % mat to hold our partial correlation data
    partialCorrMats = nan(numParcels,numParcels,numStates);
    for i = 1:numStates
        % grab timestamps, since in pairs, do 2*i + 0 or 1
        currOnTimestamps = states.(stateFields{2*i-1});
        currOffTimestamps = states.(stateFields{2*i});
        % get frames that fall within the bounds of the current state
        stateInds = boolean(zeros(numSamples,1));
        for ii = 1:length(currOffTimestamps)
            stateInds(mesoTimestamps>= currOnTimestamps(ii) & mesoTimestamps<=currOffTimestamps(ii)) = true;
        end
        % calculate partial correlation and save in standard format
        if ~isempty(currOffTimestamps)
          partialCorrMats(nanlessParcelInds,nanlessParcelInds,i) = partialcorr(nanlessParcels(:,stateInds)',meanTrace(stateInds)');
        end
    end
    saveData.partialCorrMats = partialCorrMats;
    saveData.mouse_date_time = subFolders(folder).name;
    added = added + 1;
    partialCorrData(added) = saveData;
end

for i =1:length(partialCorrData)
    figure
    for ii = 1:numStates
        subplot(numStates*100+10+ii)
        imagesc(squeeze(partialCorrData(i).partialCorrMats(:,:,ii)));
        caxis([-1 1]);
        title(partialCorrData(i).mouse_date_time)
    end
end

allPartialCorrs = nan(numParcels,numParcels,numStates,length(partialCorrData));
for i =1:length(partialCorrData)
    allPartialCorrs(:,:,:,i)=partialCorrData(i).partialCorrMats;
end


meanPartialCorrs=nanmean(allPartialCorrs,4);


stateNames = [{'Locomotion'}; {'High Face'}; {'Quiescence'}];
for i =1:numStates
        subplot(numStates*100+10+i)
        imagesc(meanPartialCorrs(:,:,i));
        caxis([-1 1]);
        title(stateNames(i))
end


desiredParcels = boolean(ones(numParcels,1));
desiredParcels(2:2:end) = false;
desiredParcels(isnan(meanPartialCorrs(1,:,1))) = false;

load('parcells_updated121519.mat')

figure
imagesc(meanPartialCorrs(desiredParcels,desiredParcels,1)-meanPartialCorrs(desiredParcels,desiredParcels,2));
caxis([-.5 .5]);
title([stateNames(1) ' - ' stateNames(2)]);
colorbar()
set(gca, 'YTick',1:sum(desiredParcels), 'YTickLabel',parcells_new.description(desiredParcels));
set(gca, 'XTick',1:sum(desiredParcels), 'XTickLabel',parcells_new.description(desiredParcels),'XTickLabelRotation',30);


figure
imagesc(meanPartialCorrs(desiredParcels,desiredParcels,1)-meanPartialCorrs(desiredParcels,desiredParcels,3));
caxis([-.5 .5]);
title([stateNames(1) ' - ' stateNames(3)]);
colorbar()
set(gca, 'YTick',1:sum(desiredParcels), 'YTickLabel',parcells_new.description(desiredParcels));
set(gca, 'XTick',1:sum(desiredParcels), 'XTickLabel',parcells_new.description(desiredParcels),'XTickLabelRotation',30);

figure
imagesc(meanPartialCorrs(desiredParcels,desiredParcels,2)-meanPartialCorrs(desiredParcels,desiredParcels,3));
caxis([-.5 .5]);
title([stateNames(2) ' - ' stateNames(3)]);
colorbar()
set(gca, 'YTick',1:sum(desiredParcels), 'YTickLabel',parcells_new.description(desiredParcels));
set(gca, 'XTick',1:sum(desiredParcels), 'XTickLabel',parcells_new.description(desiredParcels),'XTickLabelRotation',30);




imagesc(squeeze(meanPartialCorrs);