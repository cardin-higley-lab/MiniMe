function extractOddballData(oddballStimuliData, mergeOn, mergeOff)


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


