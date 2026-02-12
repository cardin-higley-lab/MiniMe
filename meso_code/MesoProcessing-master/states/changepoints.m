function [OnTStamp ,OffTStamp ] = changepoints(data, zThres,timestamps,inSampleRate,smoothWin)
%% identify change points (based on z-score, input data should be baseline normalized z-score) 
%inputs 
% data: baseline normalized z-score (pupil, face etc)
%ZThres: Z-score threshold used for detection 
%timestamps:data time in seconds 
%inSampleRate: sampling rate of data
%smoothWin: window length (in seconds) for smoothing, smoothing helps prevent triggering of Z-threshold detection by fast (noisy) fluctuations 
%b1DoPlot: whether you want to to plot the detected on/off times as 0 or 1
%b1DorPlotDuration: duration of data you want to plot
%flpPlot: if selected, flips the sign of data during plotting 
% outputs
%h1=figure; 
%OnTStamp: onset time of state
%OffTStamp: offset time of state 
%Sweyta Lohani, 2020
%% some additional parameters 
changeDur=0.5; %minimum duration of the state change in seconds 
timeBetween=0.2; % minimum time between off and the next on in seconds

%%extract time points that are above the z-threshold 
smoothDataNorm=smooth(data,smoothWin*inSampleRate)'; 
AboveThresSig          =smoothDataNorm > zThres;  %find data points greather than threshold 
tmp        =[0,AboveThresSig]; %add 0 as the first data point;  
tmp(end)  = 0; %set the last point to 0
OnIdx    = find([diff(tmp) == 1]); %on indices
OffIdx   = find([diff(tmp) == -1]); %off indices 

%if the time between off and the next on is too short, merge them 
OnIdx_tmp=OnIdx(2:end); 
OffIdx_tmp=OffIdx(1:end-1); 
idx=find((OnIdx_tmp-OffIdx_tmp)<(timeBetween*inSampleRate));
OnIdx_int=OnIdx;
OnIdx_int(idx+1)=[];
OffIdx_int=zeros(1,length(OnIdx_int)); 
for tt=1:length(OnIdx_int)-1
    indices=find(OffIdx>=OnIdx_int(tt)&OffIdx<OnIdx_int(tt+1));
    OffIdx_int(tt)=OffIdx(indices(end)); 
end 
%fill the last value
OffIdx_int(tt+1)=OffIdx(end); 

%Removes periods of changes if they are too short
wlen=round(changeDur*inSampleRate); 
RemIdx = find((OffIdx_int - OnIdx_int) <= wlen);
OnIdx_int(RemIdx) = []; 
OffIdx_int(RemIdx) = [];

%on/off timestamps 
OnTStamp   = timestamps(OnIdx_int);
OffTStamp  = timestamps(OffIdx_int);
end 