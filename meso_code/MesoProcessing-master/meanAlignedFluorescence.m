function [preEventF,postEventF] = meanAlignedFluorescence(parcelData,mesoBlueTimestamps,eventTimestamps,preEventWind,postEventWind,fs)
%% [preEventF,postEventF] = meanAlignedFluorescence(dFoF.blue,spike2_data.mesoBlueTimestamps,spike2_data.waterOnTimestamps,1,1,10)


[numParcels,~] = size(parcelData);
preEventF = nan([numParcels,length(eventTimestamps)]);
postEventF = nan(size(preEventF));

for i = 1:length(eventTimestamps)
    afterInd = find(mesoBlueTimestamps>eventTimestamps,1);
    postEventF(:,i) = mean(parcelData(:,afterInd:afterInd+postEventWind*fs));
    preEventF(:,i) = mean(parcelData(:,afterInd-1-preEventWind*fs:afterInd-1));
end
end



% mPreEventF = mean(preEventF,2);
% stdPreEventF = std(preEventF,[],2);
% 
% mPostEventF = mean(postEventF,2);
% stdPostEventF = std(postEventF,[],2);