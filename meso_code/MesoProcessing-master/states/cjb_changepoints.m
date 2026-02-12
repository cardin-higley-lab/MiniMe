function [OnTStamp ,OffTStamp] = cjb_changepoints(data,timestamps,zThresh,smoothWin)
%% identify change points (based on z-score)

%inputs 
%data: unnormalized data(pupil, face etc) (was previously normalized)
%timestamps: data time in seconds 
%smoothWin: window length (in seconds) for smoothing, smoothing helps prevent triggering of Z-threshold detection by fast (noisy) fluctuations 

%outputs
%OnTStamp: onset time of state 
%OffTStamp: offset time of state 

%Sweyta Lohani, 2020, Modified CJB 2021, modified by RLO 2021

%% some additional parameters  
timeBetween=1; % minimum time between off and the next on in seconds

if ~isrow(data) 
    data = data';
end

%% extract time points that are above the z-threshold 

AboveThres = data > zThresh;  %find data points greather than threshold 
tmp = [0,AboveThres]; %add 0 as the first data point, therefore if starts above threshold, can include;  
tmp(end) = 0; %set the last point to 0

OnIdx  = find(diff(tmp) == 1)-1; % returns indeces where it goes from 0 to 1 (below to above threshold)
if isempty(OnIdx)
    OnTStamp   = [];
    OffTStamp  = [];
    return;
end

if OnIdx(1) == 0
    OnIdx(1) = 1;
end
OffIdx = find(diff(tmp) == -1); % returns ineces where it falls below threshold, 1 to 0

%if the time between off and the next on is too short, merge them 
[OnIdx_int, OffIdx_int] = mergeTimestamps(OnIdx, OffIdx, timeBetween);


outsideRange = find(OffIdx_int > length(timestamps));
if ~isempty(outsideRange)
    OffIdx_int(outsideRange) = [];
    OnIdx_int = OnIdx_int(1:length(OffIdx_int));
end


%on/off timestamps 
OnTStamp   = timestamps(OnIdx_int)';
OffTStamp  = timestamps(OffIdx_int)';


end 