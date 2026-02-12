[preEventF,postEventF] = meanAlignedFluorescence(dFoF.blue,spike2_data.mesoBlueTimestamps,spike2_data.waterOnTimestamps,1,1,10);

mPreEventF = mean(preEventF,2);
stdPreEventF = std(preEventF,[],2);

mPostEventF = mean(postEventF,2);
stdPostEventF = std(postEventF,[],2);

numParcels = 56;
pVals = zeros([numParcels,1]);
for i = 1:length(numParcels)
    [h,p] = ttest(preEventF(i,:),postEventF(i,:));
    pVals(i) = p;
end

[FDR,Q] = mafdr(pVals); 