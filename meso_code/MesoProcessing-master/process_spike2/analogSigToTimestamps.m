function  [onTimestamps, offTimestamps] = analogSigToTimestamps(inSig,inTime)
%%analogSigToTimestamps transform 2 state analog signal to timestamps 
%   Employs thresholds halfway between min and max values to convert an analog signal into timestamps.
%   inSig: analog signal states are extracted from
%   inTime: time vector used to create timestamps for state transitions


if mean(inSig,'omitnan') < 0
    inSig = -1*inSig;
end

stateVect = inSig>nanmin(inSig)+0.5*(nanmax(inSig)-nanmin(inSig)); % nanmean + thresh to fix weird negative voltage visual triggers

stateChange = diff(stateVect); % transform into state changes

onTimestamps = inTime(find(stateChange==1)+1); % get time of state changes
offTimestamps = inTime(find(stateChange==-1)+1);  % get time of state changes
if ~isempty(onTimestamps) && ~isempty(offTimestamps)
    if onTimestamps(1) > offTimestamps(1) % if assigned states "wrong", flip on and off timestamps
        tempTimestamps = onTimestamps;
        onTimestamps = offTimestamps;
        offTimestamps = tempTimestamps;
    end
end

if length(onTimestamps)-length(offTimestamps)==1 % program stopped mid state, so delete last rising edge
    onTimestamps(end)=[];
end
end

