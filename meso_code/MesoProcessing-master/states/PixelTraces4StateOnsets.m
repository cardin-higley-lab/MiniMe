function pixelwise_stateonset = PixelTraces4StateOnsets(samplingrate, preSeconds, postSeconds, states, spike2, normed_sig)


    mesoTimestamps = spike2.blueOnTimestamps; 
    numSamples = size(normed_sig, 2);
    maxLength= min(length(mesoTimestamps), numSamples);  % make sure theyre all the same length
    mesoTimestamps = mesoTimestamps(1:maxLength);
    normed_sig = normed_sig(:, 1:maxLength);


    preSamples = round(samplingrate*preSeconds); % convert seconds to num samples (frames)
    postSamples = round(samplingrate*postSeconds);
    numFrames = preSamples + 1 + postSamples;
    numPixels = size(normed_sig, 1);

    stateFields = fieldnames(states);
    numStates = length(stateFields)/2;

    pixelwise_stateonset = zeros(numPixels, numFrames, numStates, 'single');
    for state = 1:numStates
  
        state_name = stateFields{(state*2)-1};
        stateOn = states.(state_name);
        numBouts = length(stateOn);

        if state_name == "locoOn"
            stateOn = stateOn - 3;
        elseif state_name == "sitOn"
            stateOn = stateOn - 10;
        end

        %pre-allocations
        meanPixels = nan(numBouts, numPixels, numFrames); % for events, parcels, and total frames
        numValidEvents = 0; % counts number of valid events, used to calculate standard error of the mean (SEM)

        % iterate through state events to find associated meso frame data
        for bout = 1:numBouts
  
            ind = find(mesoTimestamps>=stateOn(bout),1); % get first meso frame on or after event time
            if ~isempty(ind) > 0
                if ind > 1 % since we will subtract 1 from ind, make sure it will not become 0
                    possibleInds = [ind-1 ind]; % since found first after, previous ind could be closer
                    [~,bestIndInd] = min([mesoTimestamps(ind-1)-stateOn(bout) mesoTimestamps(ind)-stateOn(bout)]); % check which ind is closer to the event
                    bestInd = possibleInds(bestIndInd); % set the best ind
                else
                    bestInd = ind; % ind is 1
                end

                if bestInd-preSamples>0 && bestInd+postSamples<=numSamples % check to see if event occured in a manner we can align with meso data

                    meanPixels(bout,:,:) = normed_sig(:,bestInd-preSamples:bestInd+postSamples); % grab the trace
                    numValidEvents = numValidEvents + 1; % add another valid event
                end
            end

        end

       
        pixelMean = squeeze(mean(meanPixels,1, 'omitnan')); % take eventwise mean, ignoring nans (nan = couldnt align meso). Squeeze out now unecessary dimension
        pixelwise_stateonset(:,:,state) = pixelMean;
    end
 

end 