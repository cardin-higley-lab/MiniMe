function  [onTimestamps, offTimestamps] = analog_sig_to_timestamps(inSig,inTime)
%%analog_sig_to_timestamps transform 2 state analog signal to timestamps 
%   Employs kmeans clustering or simple thresholding, dependent on signal std, to convert an analog signal into timestamps.
%   inSig: analog signal states are extracted from
%   inTime: time vector used to create timestamps for state transitions

%stateVect = (inSig-nanmin(inSig))>0.5;
%if nanstd(inSig) > 0.5 % signal is substantially variable (appropriate for kmeans)
%    stateVect = kmeans(inSig,2);
%else % signal is sparse (not appropriate for kmeans clustering)
    inSig = inSig - nanmin(inSig);
    stateVect = inSig>0.5*nanmax(inSig);
%end

stateChange = diff(stateVect); % transform into state changes
onTimestamps = inTime(stateChange==1); % get time of state changes
offTimestamps = inTime(stateChange==-1);  % get time of state changes
if ~isempty(onTimestamps) && ~isempty(offTimestamps)
    if onTimestamps(1) > offTimestamps(1) % if kmeans assigned states "wrong", flip on and off timestamps
        tempTimestamps = onTimestamps;
        onTimestamps = offTimestamps;
        offTimestamps = tempTimestamps;
    end
end

if length(onTimestamps)-length(offTimestamps)==1 % program stopped mid state, so delete last rising edge
    onTimestamps(end)=[];
end
end

