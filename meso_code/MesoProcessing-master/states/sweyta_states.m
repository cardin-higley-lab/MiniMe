function sweyta = sweyta_states(spike2_data, face_data, hands)

params.fspupilcam=10; %pupil sampling rate
params.fsspike2=5000;% spike2 sampling rate
params.TimeSinceLocOn=3;%for locomotion state, minimum time since locomotion onset
params.TimeBeforeLocOff=3;%for locomotion state, minimum time before locomotion offset
params.TimeSinceSitOn=10;%for quiescence state, minimum time since quiescence onset
params.TimeBeforeSitOff=10;%for quiescence state, minimum time before quiescence offset
params.TimeSinceEvent=10;%for any state, minimum time since any event onset/offset
params.minRunDuration=2;% minimum run duration during locomotion state, previous 5s now 2s
params.minArousalDuration=2; %minimum face/pupil arousal state (high or low arousal),previous 5s now 2s
params.minSitDuration=2;%minimum sit duration during quiescnece state,previous 5s now 2s


face_Norm=normalize(face_data);
wheel_speed=spike2_data.wheelSpeed./100;
% get pupil wheel and imaging times
wheel_time = (1:length(wheel_speed))/params.fsspike2;
imaging_time=spike2_data.blueOnTimestamps;
pupil_time = spike2_data.pupilFrameOnTimestamps (1:length(face_Norm));
%% get locomotion on/off and quiescence on off times
%locomotion periods should be at least some criterion s long with some criterion sec since locomotion onset, some criterion sec before locomotion offset, excluding any events (airpuff/stim)
if spike2_data.wheelOn(1) > spike2_data.wheelOff(1) % account for animal running when recording started
    wheelOn = [0; spike2_data.wheelOn];
else 
    wheelOn = spike2_data.wheelOn;
end

if spike2_data.wheelOff(end) < spike2_data.wheelOn(end) % account for animal running when recording started
    wheelOff = [spike2_data.wheelOff; pupil_time(end)];
else 
    wheelOff = spike2_data.wheelOff;
end

%find wheel on and wheel off times during imaging period only
minRunDur=params.minRunDuration+params.TimeSinceLocOn+params.TimeBeforeLocOff; %minimum actual locomotion duration including time since locomotion onset, time before locomotion offset and the minimum time period for data analysis
idx=wheelOn<(imaging_time(end)) & wheelOff>(imaging_time(1));
wheelOn_t1=wheelOn(idx);
wheelOff_t1=wheelOff(idx);

for whe=1:length(wheelOn_t1)
    if wheelOn_t1(whe)<imaging_time(end) && wheelOff_t1(whe)>imaging_time(end)
        wheelOff_t1(whe)=imaging_time(end)-1;%if locomotion starts before end of imaging but continues after, only extract state until imaging time end minus a second
    end
    
    if wheelOn_t1(whe)<imaging_time(1) && wheelOff_t1(whe)>imaging_time(1)
        wheelOn_t1(whe)=imaging_time(1)+1;%if locomotion starts before start of imaging but continues after, only extract state from imaging time start plus a second
    end
end

wheelOn_t1=(wheelOn_t1(:))'; wheelOff_t1=(wheelOff_t1(:))';


wheelOn_int=wheelOn_t1;
wheelOff_int=wheelOff_t1;

%makes sure the state is at least as long as the minimum run duration
idx1=find((wheelOff_int-wheelOn_int)>=(minRunDur));
wheelOn_int1=wheelOn_int(idx1);
wheelOff_int1=wheelOff_int(idx1);

%finalize the times to get sustained state only
wheelOn_final=wheelOn_int1+params.TimeSinceLocOn;
wheelOff_final=wheelOff_int1-params.TimeBeforeLocOff;

sweyta.locoOn = wheelOn_final;
sweyta.locoOff = wheelOff_final;

%% queiscence should be at least some criterion s long with some criterion s since locomotion offset
%and some criterion s before subsequent locomotion onset,excluding any events (airpuff/stim)
sitOn=[0;wheelOff]; %use 0 as the first sit on time;
sitOff=[wheelOn;imaging_time(end)];%use wheelOn times as sit off times;


%find sit on and sit off times during imaging period only
minSitDur=params.minSitDuration+params.TimeSinceSitOn+params.TimeBeforeSitOff; %actual minimum sit duration accouting for the onset time, offset time and minimum duration of the sustained quiescence epoch  used for analysis

idx=sitOn<(imaging_time(end)) & sitOff>(imaging_time(1));
sitOn_t1=sitOn(idx);
sitOff_t1=sitOff(idx);

for whe=1:length(sitOn_t1)
    if sitOn_t1(whe)<imaging_time(end) && sitOff_t1(whe)>imaging_time(end)
        sitOff_t1(whe)=imaging_time(end)-1;%if queiscence starts before end of imaging but continues after, only extract state until imaging time end minus a second
    end
    
    if sitOn_t1(whe)<imaging_time(1) && sitOff_t1(whe)>imaging_time(1)
        sitOn_t1(whe)=imaging_time(1)+1;%if queiscence starts before start of imaging but continues after, only extract state from imaging time start plus a second
    end
end

wheel_speed=wheel_speed(:);
wheel_time=wheel_time(:);
%remove any quiescence period where mouse's speed is above 0.03 m/s
% because sometimes these epcohs,especially short ones, get missed by locomotion changepoint algorithm
speedThres=0.03;
tmpOn=sitOn_t1; tmpOff=sitOff_t1;
for rr=1:length(sitOn_t1)
    highSpeedIdx=find(abs(wheel_speed)>speedThres & wheel_time>sitOn_t1(rr) & wheel_time<sitOff_t1(rr));
    if ~isempty(highSpeedIdx)
        firstIdx=wheel_time(highSpeedIdx(1))-0.1; lastIdx=wheel_time(highSpeedIdx(end))+0.1;
        tmpOff(rr)=firstIdx; tmpOff(end+1)=sitOff_t1(rr); tmpOn(end+1)=lastIdx;
    end
end
tmpOn=sort(tmpOn,'ascend'); tmpOff=sort(tmpOff,'ascend');

sitOn_t1=(tmpOn(:))'; sitOff_t1=(tmpOff(:))';


sitOn_int=sitOn_t1;
sitOff_int=sitOff_t1;

%makes sure state is at least as long as the minimum sit duration
idx1=find((sitOff_int-sitOn_int)>=(minSitDur));
sitOn_int1=sitOn_int(idx1);
sitOff_int1=sitOff_int(idx1);

%finalize the times to get sustained state only
sitOn_final=sitOn_int1+params.TimeSinceSitOn;
sitOff_final=sitOff_int1-params.TimeBeforeSitOff;

sweyta.sitOn = sitOn_final;
sweyta.sitOff = sitOff_final;



%% do change point detection on face to get face high/low movement times during sustained quiescence state
%get Z-thresholds based on face data during quiescence, when mouse isn't moving and when aripuffs are not given


% gets time points within sitting bouts 
pupilTime_Idx=cell(1,length(sitOn_int));
for st=1:length(sitOn_int)
    pupilTime_Idx{st}=find(pupil_time>sitOn_int(st) & pupil_time <sitOff_int(st));
end
pupilTime_quiescence=cell2mat(pupilTime_Idx');

% gets normalized face trace during quiescent outs (if not sit, then takes
% all of normalized face trace)
if ~isempty(pupilTime_quiescence)
    face_quiescence=face_Norm(pupilTime_quiescence);
else
    face_quiescence=face_Norm;
end
zthres_High=quantile(face_quiescence,0.60);
zthres_Low=quantile(face_quiescence,0.40);

%get on and off timestamps for high and low face movment
% [Face_HighArousal_OnTStamp,Face_HighArousal_OffTStamp] = changepoints(face_Norm', zthres_High,pupil_time,params.fspupilcam, 1);
% [Face_LowArousal_OnTStamp,Face_LowArousal_OffTStamp] = changepoints(-face_Norm', -zthres_Low,pupil_time,params.fspupilcam, 1);

% with Clayton's function
[Face_HighArousal_OnTStamp,Face_HighArousal_OffTStamp] = cjb_changepoints(face_Norm, pupil_time, zthres_High,1);
[Face_LowArousal_OnTStamp,Face_LowArousal_OffTStamp] = cjb_changepoints(-face_Norm, pupil_time, -zthres_Low,1); 



%remove outliers (due to grooming for example) in high face state data by removing states where the state values are greater than 2 standard deviation from the whole session average
znormedface=normalize(face_Norm);

maxzdata=nan(1,length(Face_HighArousal_OnTStamp));
for st=1:length(Face_HighArousal_OnTStamp)
    OnIndx=find(pupil_time==Face_HighArousal_OnTStamp(st));
    OffIndx=find(pupil_time==Face_HighArousal_OffTStamp(st));
    zdata=znormedface(OnIndx:OffIndx);
    maxzdata(st)=nanmax(zdata);
end

Face_HighArousal_OnTStamp=Face_HighArousal_OnTStamp(maxzdata<4);
Face_HighArousal_OffTStamp=Face_HighArousal_OffTStamp(maxzdata<4);

% get  face on times if both on and off times occur during
% sustained queiscence states identified in the previous step. If on/off times occur beyond the quiescence state, modify the on/off times to be during the quiescence only
s1=sitOn_final; 
s2=sitOff_final;
on_final=cell(1,length(Face_HighArousal_OnTStamp)); 
off_final=cell(1,length(Face_HighArousal_OnTStamp));
for rj=1:length(Face_HighArousal_OnTStamp)
    a1=Face_HighArousal_OnTStamp(rj); a2=Face_HighArousal_OffTStamp(rj);
    A=a1-s1; B=a2-s2;
    on_f=nan(1,length(sitOn_final)); off_f=nan(1,length(sitOff_final));
    for rt=1:length(sitOn_final)
        if A(rt)>0,Ac=A(rt); else Ac=0; end
        if B(rt)<0,Bc=abs(B(rt)); else Bc=0; end
        on= s1(rt)+Ac;
        off=s2(rt)-Bc;
        if (off-on)>0,on_f(rt)=on; off_f(rt)=off; end
    end
    on_final{rj}=on_f(find(~isnan(on_f))); off_final{rj}=off_f(find(~isnan(on_f)));
end
Face_HighArousalOn_int1=cell2mat(on_final); Face_HighArousalOff_int1=cell2mat(off_final);

on_final=cell(1,length(Face_LowArousal_OnTStamp)); off_final=cell(1,length(Face_LowArousal_OnTStamp));
for rj=1:length(Face_LowArousal_OnTStamp)
    a1=Face_LowArousal_OnTStamp(rj); a2=Face_LowArousal_OffTStamp(rj);
    A=a1-s1; B=a2-s2;
    on_f=nan(1,length(sitOn_final)); off_f=nan(1,length(sitOff_final));
    for rt=1:length(sitOn_final)
        if A(rt)>0,Ac=A(rt); else Ac=0; end
        if B(rt)<0,Bc=abs(B(rt)); else Bc=0; end
        on= s1(rt)+Ac;
        off=s2(rt)-Bc;
        if (off-on)>0,on_f(rt)=on; off_f(rt)=off; end
    end
    on_final{rj}=on_f(find(~isnan(on_f))); off_final{rj}=off_f(find(~isnan(on_f)));
end
Face_LowArousalOn_int1=cell2mat(on_final); Face_LowArousalOff_int1=cell2mat(off_final);

% determine that face high/low arousal/movement time are at least minimum criterion seconds long
idx3=find((Face_HighArousalOff_int1-Face_HighArousalOn_int1)>=params.minArousalDuration);
Face_HighArousal_On_final=Face_HighArousalOn_int1(idx3); Face_HighArousal_Off_final=Face_HighArousalOff_int1(idx3);

idx3=find((Face_LowArousalOff_int1-Face_LowArousalOn_int1)>=params.minArousalDuration);
Face_LowArousal_On_final=Face_LowArousalOn_int1(idx3); Face_LowArousal_Off_final=Face_LowArousalOff_int1(idx3);


sweyta.faceHighSitOn = Face_HighArousal_On_final; 
sweyta.faceHighSitOff = Face_HighArousal_Off_final;


sweyta.faceLowSitOn = Face_LowArousal_On_final; 
sweyta.faceLowSitOff = Face_LowArousal_Off_final;


%% plot on and off times for sanity check
%plot locomotion state on off times
h=figure; %initialize a figure that will contain subplots of behavior with on/off times markers for identified states
set(0,'CurrentFigure',h);ax1=subplot(4,1,1);plot(wheel_time,wheel_speed);ylimits=ylim; hold on; title('LocomotionState');
for tt=1:length(wheelOn_final)
    plot([wheelOn_final(tt),wheelOn_final(tt)], ylimits, 'g');
    plot([wheelOff_final(tt),wheelOff_final(tt)], ylimits, 'r');
end

plot([imaging_time(1),imaging_time(1)], ylimits, 'm');
plot([imaging_time(end),imaging_time(end)], ylimits, 'm');

%plot quiescence state on off times
set(0,'CurrentFigure',h);ax2=subplot(4,1,2);plot(wheel_time,wheel_speed);ylimits=ylim; hold on; title('QuiescenceState');
for tt=1:length(sitOn_final)
    plot([sitOn_final(tt),sitOn_final(tt)], ylimits, 'g');
    plot([sitOff_final(tt),sitOff_final(tt)], ylimits, 'r');
end

plot([imaging_time(1),imaging_time(1)], ylimits, 'm');
plot([imaging_time(end),imaging_time(end)], ylimits, 'm');

%plot face state on off times
set(0,'CurrentFigure',h);ax5=subplot(4,1,3);plot(pupil_time,face_Norm);ylimits=ylim; hold on; title('FaceHighState');
for tt=1:length(Face_HighArousal_On_final)
    plot([Face_HighArousal_On_final(tt),Face_HighArousal_On_final(tt)], ylimits, 'g');
    plot([Face_HighArousal_Off_final(tt),Face_HighArousal_Off_final(tt)], ylimits, 'r');
end

set(0,'CurrentFigure',h);ax6=subplot(4,1,4);plot(pupil_time,face_Norm);ylimits=ylim; hold on; title('FaceLowState');
for tt=1:length(Face_LowArousal_On_final)
    plot([Face_LowArousal_On_final(tt),Face_LowArousal_On_final(tt)], ylimits, 'g');
    plot([Face_LowArousal_Off_final(tt),Face_LowArousal_Off_final(tt)], ylimits, 'r');
end

linkaxes([ax1, ax2,ax5,ax6],'x');

%% adding grooming 

minStateDur = 2; %making all states need to be at least 5 sec to count 
TimeSinceSitOn=10;%for quiescence state, minimum time since quiescence onset
TimeBeforeSitOff=10;%for quiescence state, minimum time before quiescence offset


cutoff_freq = 0.1;
sampling_freq = round(1/median(diff(spike2_data.pupilFrameOnTimestamps)));
[b,a] = butter(2,cutoff_freq/(sampling_freq/2));


if std(hands) > 200
    hands = filtfilt(b,a,double(hands));
    hands = zscore(hands);
    hands_threshold = max(hands)/2;

    quant_thresh=quantile(hands, 0.6);
    
    [groomHighOn, groomHighOff] = cjb_changepoints(hands,spike2_data.pupilFrameOnTimestamps,quant_thresh,1); 
    groomHighOn_int=squeeze(groomHighOn(groomHighOff-groomHighOn>minStateDur));
    groomHighOff_int=squeeze(groomHighOff(groomHighOff-groomHighOn>minStateDur));
    
    
    % remove any grooming period where z-score is not above 1
    groomHighOn_int2 = zeros(1, length(groomHighOn_int));
    groomHighOff_int2 = zeros(1, length(groomHighOff_int));
    
    
    for groom_bout = 1:length(groomHighOn_int)
        ind1 = find(spike2_data.pupilFrameOnTimestamps == groomHighOn_int(groom_bout));
        ind2 = find(spike2_data.pupilFrameOnTimestamps == groomHighOff_int(groom_bout));
        groom_max = max(hands(ind1:ind2-1));
        if groom_max > hands_threshold
            groomHighOn_int2(groom_bout) = groomHighOn_int(groom_bout); 
            groomHighOff_int2(groom_bout) = groomHighOff_int(groom_bout);
    
        end
        if groom_max <= hands_threshold
            groomHighOn_int2(groom_bout) = NaN;
            groomHighOff_int2(groom_bout) = NaN;
    
        end
    end
    groomHighOn_final = groomHighOn_int2(~isnan(groomHighOn_int2));
    groomHighOff_final = groomHighOff_int2(~isnan(groomHighOff_int2));
    
    
    sweyta.groomHighOn = groomHighOn_final ; 
    sweyta.groomHighOff = groomHighOff_final ;

elseif std(hands) <= 200
    sweyta.groomHighOn = [];
    sweyta.groomHighOff = [];

end 




end