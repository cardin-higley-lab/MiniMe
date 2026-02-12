function [OnTStamp ,OffTStamp] = cjb_changepoints_bw(data,timestamps,zThresh_upper,zThresh_lower, smoothWin)
%% identify change points (based on z-score)

%inputs 
%data: unnormalized data(pupil, face etc) (was previously normalized)
%timestamps: data time in seconds 
%quantileThresh: quantile to normalize data [0 1] 
%smoothWin: window length (in seconds) for smoothing, smoothing helps prevent triggering of Z-threshold detection by fast (noisy) fluctuations 

%outputs
%OnTStamp: onset time of state 
%OffTStamp: offset time of state 

%Sweyta Lohani, 2020, Modified CJB 2021, modified by RLO 2021

%% some additional parameters 
changeDur=0.5; %minimum duration of the state in seconds 
timeBetween=0.2; % minimum time between off and the next on in seconds
inSampleRate = 1/median(diff(timestamps)); % finds the average sampling rate 

% % normalize data (& and smooth for theshold detection)
% Ms = sort(data(~isnan(data)),'ascend');% Sort asending along time dimension
% F0 = Ms(1:ceil(length(Ms)*0.1)); % lower 10% of the values
% MeanF0=mean(F0);
% normData=(data-MeanF0)./std(F0); % zscoring.. uses a diff mean than base zscore fx would

normData_sm = smooth(data,smoothWin*inSampleRate)'; 

%% extract time points that are above the z-threshold 
withinThres = zThresh_upper > normData_sm & normData_sm > zThresh_lower;  %find data points greather than threshold 
tmp = [0,withinThres]; %add 0 as the first data point, therefore if starts above threshold, can include;  
tmp(end) = 0; %set the last point to 0
OnIdx  = find(diff(tmp) == 1)+1; % returns indeces where it goes from 0 to 1 (below to above threshold)
OffIdx = find(diff(tmp) == -1)+1; % returns ineces where it falls below threshold, 1 to 0

%if the time between off and the next on is too short, merge them 
OnIdx_tmp=OnIdx(2:end); % no previous off to potentiall merge with
OffIdx_tmp=OffIdx(1:end-1); % no next on  to potentially merge with
idx=find((OnIdx_tmp-OffIdx_tmp)<(timeBetween*inSampleRate)); %finds where state dur is too short
OnIdx_int=OnIdx;
OnIdx_int(idx+1)=[]; %add 1 back bc removed to find 'idx'; removing these, therefore merging the periods
OffIdx_int=zeros(1,length(OnIdx_int)); 
for tt=1:length(OnIdx_int)-1
    indices=find(OffIdx>=OnIdx_int(tt)&OffIdx<OnIdx_int(tt+1));
    OffIdx_int(tt)=OffIdx(indices(end)); 
end 
OffIdx_int(tt+1)=OffIdx(end); %fill the last value

%Removes periods of changes if they are too short
window_length=round(changeDur*inSampleRate); 
RemIdx = find((OffIdx_int - OnIdx_int) <= window_length);
OnIdx_int(RemIdx) = []; 
OffIdx_int(RemIdx) = [];

outsideRange = find(OffIdx_int > length(timestamps));
if ~isempty(outsideRange)
    OffIdx_int(outsideRange) = [];
    OnIdx_int = OnIdx_int(1:length(OffIdx_int));
end

%on/off timestamps 
OnTStamp   = timestamps(OnIdx_int);
OffTStamp  = timestamps(OffIdx_int);
end 