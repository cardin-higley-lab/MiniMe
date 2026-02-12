function states = stateTimestamps(spike2_data,pc, hands, varargin)
%STATETIMESTAMPS identifies states from wheel and facial motion
%   returns timetamps for high face no run, run (almost always accompanied by high face), and quiescence (face low
%   and no running). 
% Variable arguments allow for changing default parameters:
%   ex: stateTimestamps(spike2_data,pc,'minStateDuration',3)
%   this will change the minStateDur to 3 seconds from 5 (default)
% Assumes: 
% 1) facemap .mat file ends in proc.mat
% 2) only 1 motion SVD region was selected
% 

minStateDur = 5; %making all states need to be at least 5 sec to count 
TimeSinceLocOn=3;%for locomotion state, minimum time since locomotion onset
TimeBeforeLocOff=3;%for locomotion state, minimum time before locomotion offset
TimeSinceSitOn=10;%for quiescence state, minimum time since quiescence onset
TimeBeforeSitOff=10;%for quiescence state, minimum time before quiescence offset
fsspike2=5000;

%% change default values
if any(strcmp('minStateDuration',varargin))
    disp('minStateDuration')
    minStateDuration = varargin{find(strcmp('minStateDuration',varargin))+1};
end
if any(strcmp('timeSinceLocOn',varargin))
	timeSinceLocOn = varargin{find(strcmp('timeSinceLocOn',varargin))+1};
end
if any(strcmp('timeBeforeLocOff',varargin))
    timeBeforeLocOff = varargin{find(strcmp('timeBeforeLocOff',varargin))+1};
end
if any(strcmp('timeSinceSitOn',varargin))
    timeSinceSitOn = varargin{find(strcmp('timeSinceSitOn',varargin))+1};
end
if any(strcmp('timeBeforeSitOff',varargin))
    timeBeforeSitOff = varargin{find(strcmp('timeBeforeSitOff',varargin))+1};
end
if any(strcmp('fsspike2',varargin))
	fsspike2 = varargin{find(strcmp('fsspike2',varargin))+1};
end
if any(strcmp('raw',varargin))
    minStateDuration = 0;
    timeSinceLocOn = 0;
    timeBeforeLocOff = 0;
    timeSinceSitOn = 0;
    timeBeforeSitOff = 0;
end


%% Locomotion 

if spike2_data.wheelOn(1) > spike2_data.wheelOff(1) % account for animal running when recording started
    wheelOff = spike2_data.wheelOff(2:end);
else 
    wheelOff = spike2_data.wheelOff;
end
wheelOn = spike2_data.wheelOn(1:length(wheelOff)); % account for animal running when recording stopped


blueMesoTimestamps = spike2_data.blueOnTimestamps;

% min loco duration accouting for the onset/offset times and min state period for data analysis
minRunDur=minStateDur+TimeSinceLocOn+TimeBeforeLocOff;

firstLocoOn = wheelOn-minRunDur>=blueMesoTimestamps(1); %dont necessarily need to the -minRunDur
lastLocoOff= wheelOff+minRunDur<=blueMesoTimestamps(end); %dont necessarily need to the +minRunDur
locoBoutTime = wheelOff -wheelOn > minRunDur; 

wheelOn_final = wheelOn(firstLocoOn & lastLocoOff & locoBoutTime);
wheelOff_final = wheelOff(firstLocoOn & lastLocoOff & locoBoutTime);

% add to states data
states.locoOn = wheelOn_final' + TimeSinceLocOn;
states.locoOff = wheelOff_final' - TimeBeforeLocOff;

%% Sit (Quiescience)

wheelSpeed = spike2_data.wheelSpeed;
wheelTime = (1:length(wheelSpeed))/fsspike2;
sitOn=[0 reshape(wheelOff,1,[])]; %use 0 as the first sit on time;
sitOff=[reshape(wheelOn,1,[]) wheelTime(end)];%use wheelOn times as sit off times;

% min sit duration accouting for the onset/offset times and min state period for data analysis
minSitDur=minStateDur+TimeSinceSitOn+TimeBeforeSitOff; 

firstSitOn = sitOn-minSitDur>=blueMesoTimestamps(1);
lastSitOff= sitOff+minSitDur<=blueMesoTimestamps(end);
sitBoutTime = sitOff - sitOn > minSitDur; 

%find sit on and sit off times during imaging period only
sitOn_final = sitOn(firstSitOn & lastSitOff & sitBoutTime);
sitOff_final = sitOff(firstSitOn & lastSitOff & sitBoutTime);

% add to states data
states.sitOn = sitOn_final + TimeSinceSitOn;
states.sitOff = sitOff_final - TimeBeforeSitOff; 

%% Grooming

cutoff_freq = 0.1;
sampling_freq = round(1/median(diff(spike2_data.pupilFrameOnTimestamps)));
[b,a] = butter(2,cutoff_freq/(sampling_freq/2));
hands = filtfilt(b,a,double(hands));

[groomHighOn, groomHighOff] = cjb_changepoints(hands,spike2_data.pupilFrameOnTimestamps,0.6,1);
groomHighOn_inter=squeeze(groomHighOn(groomHighOff-groomHighOn>minSitDur));
groomHighOff_inter=squeeze(groomHighOff(groomHighOff-groomHighOn>minSitDur));

% face high & wheel low 
[groomHighSitOn,groomHighSitOff] = timestampsContained(groomHighOn_inter,groomHighOff_inter,sitOn_final,sitOff_final);

states.groomHighSitOn = groomHighSitOn + TimeSinceSitOn; 
states.groomHighSitOff = groomHighSitOff - TimeBeforeSitOff;

%% for visualizing
figure()
maxlength = min(length(spike2_data.pupilFrameOnTimestamps), length(hands));
plot(spike2_data.pupilFrameOnTimestamps(1:maxlength), hands(1:maxlength));
xline(groomHighOn_inter, 'g')
xline(groomHighOff_inter, 'r')

figure()
maxlength = min(length(spike2_data.pupilFrameOnTimestamps), length(hands));
plot(spike2_data.pupilFrameOnTimestamps(1:maxlength), hands(1:maxlength));
xline(groomHighSitOn, 'g')
xline(groomHighSitOff, 'r')


%% HIGH FACE MAP

pc = filtfilt(b,a,double(pc));


[faceHighOn,faceHighOff] = cjb_changepoints(pc,spike2_data.pupilFrameOnTimestamps,0.6,1);
faceHighOn_inter=squeeze(faceHighOn(faceHighOff-faceHighOn>minSitDur));
faceHighOff_inter=squeeze(faceHighOff(faceHighOff-faceHighOn>minSitDur));

% face high & wheel low 
[faceHighSitOn,faceHighSitOff] = timestampsContained(faceHighOn_inter,faceHighOff_inter,sitOn_final,sitOff_final);

states.faceHighSitOn = faceHighSitOn + TimeSinceSitOn; 
states.faceHighSitOff = faceHighSitOff - TimeBeforeSitOff;


%% LOW FACE MAP

[faceLowOn,faceLowOff] = cjb_changepoints(-pc,spike2_data.pupilFrameOnTimestamps,0.4,1); 

[faceLowSitOn,faceLowSitOff] = timestampsAnd(faceLowOn,faceLowOff,sitOn_final,sitOff_final);
faceLowSitOn_inter = squeeze(faceLowSitOn(faceLowSitOff-faceLowSitOn>minSitDur));
faceLowSitOff_inter = squeeze(faceLowSitOff(faceLowSitOff-faceLowSitOn>minSitDur));

states.faceLowSitOn = faceLowSitOn_inter + TimeSinceSitOn;
states.faceLowSitOff = faceLowSitOff_inter - TimeBeforeSitOff;




end

