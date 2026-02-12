function parcelMean  = parcels4locoOnset(dff, blueOn, blueOff, wheelOn, fs, preSeconds, postSeconds)    

%     load([fileLoc, '\final_dFoF_parcels.mat']); % load parcel data
%     load([fileLoc, '\final_timestamps.mat']); % load timestamps
%     
    [numParcels,numSamples] = size(dff); % get the total number of parcels and samples
    
    % get mid point of meso frames (time point most representative of data)
    maxLength = min(length(blueOn), length(blueOff));
    mesoTimestamps = (blueOff(1:maxLength)+ blueOn(1:maxLength))/2;
    
    %preallocations:
    preSamples = round(fs*preSeconds); % convert seconds to num samples (frames)
    postSamples = round(fs*postSeconds);
    meanParcels = nan(length(wheelOn),numParcels,preSamples+postSamples+1); % for events, parcels, and total frames
    numValidEvents = 0; % counts number of valid events, used to calculate standard error of the mean (SEM)
    
    % iterate through wheelOn events to find associated meso frame data
    for i = 1:length(wheelOn) % iterate through events
        ind = find(mesoTimestamps>=wheelOn(i),1); % get first meso frame on or after event time
        if ~isempty(ind) > 0
            if ind > 1 % since we will subtract 1 from ind, make sure it will not become 0
                possibleInds = [ind-1 ind]; % since found first after, previous ind could be closer
                [~,bestIndInd] = min([mesoTimestamps(ind-1)-wheelOn(i) mesoTimestamps(ind)-wheelOn(i)]); % check which ind is closer to the event
                bestInd = possibleInds(bestIndInd); % set the best ind
            else
                bestInd = ind; % ind is 1
            end
            if bestInd-preSamples>0 && bestInd+postSamples<=numSamples % check to see if event occured in a manner we can align with meso data
                meanParcels(i,:,:) = dff(:,bestInd-preSamples:bestInd+postSamples); % grab the traces
             
                numValidEvents = numValidEvents + 1; % add another valid event
            end
        end
    end
    
    parcelMean = squeeze(nanmean(meanParcels,1)); % take eventwise mean, ignoring nans (nan = couldnt align meso). Squeeze out now unecessary dimension
   

end