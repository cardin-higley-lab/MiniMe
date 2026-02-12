function [alignedData, preData] = averageMovieToEvent(data,mesoTimestamps,eventTimestamps,preSamples,postSamples)

[numPixels,t]= size(data);

minSamples = min(t,length(mesoTimestamps));
mesoTimestamps = mesoTimestamps(1:minSamples);
data = data(:,1:minSamples);
alignedData = nan(length(eventTimestamps),numPixels,preSamples+postSamples, 'single');  
preData = nan(length(eventTimestamps),numPixels,preSamples, 'single');  


for i = 1:length(eventTimestamps)
    ind = findInSorted(mesoTimestamps,eventTimestamps(i));
    if ind-preSamples> 0 && ind+postSamples-1 <= minSamples
        alignedData(i,:,:) = data(:,ind-preSamples:ind+postSamples-1);
        preData(i,:,:) = data(:,ind-preSamples:ind-1);
    end
end