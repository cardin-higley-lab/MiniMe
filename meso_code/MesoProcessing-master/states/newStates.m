function states = newStates(spike2_data, vid_energy, state_duration, loco_buffer, sit_buffer, session)
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
% suggested time durations: state_duration=2; loco_buffer=3; sit_buffer=10;
 
   
minStateDur = state_duration; %making all states need to be at least 5 sec to count 
minGroomDur = 5; 
TimeSinceLocOn=loco_buffer;%for locomotion state, minimum time since locomotion onset
TimeBeforeLocOff=loco_buffer;%for locomotion state, minimum time before locomotion offset
TimeSinceSitOn=sit_buffer;%for quiescence state, minimum time since quiescence onset
TimeBeforeSitOff=sit_buffer;%for quiescence state, minimum time before quiescence offset
fsspike2=5000;

minSitDur=minStateDur+TimeSinceSitOn+TimeBeforeSitOff; 
minRunDur=minStateDur+TimeSinceLocOn+TimeBeforeLocOff;
max_pupil = min(length(spike2_data.pupilFrameOnTimestamps), length(vid_energy));
pupil_time = spike2_data.pupilFrameOnTimestamps(1:max_pupil);
wheelTime = (1:length(spike2_data.wheelSpeed))/fsspike2;

if isfield(spike2_data, 'blueOnTimestamps')
    imaging_time = spike2_data.blueOnTimestamps;
else
    imaging_time = spike2_data.greenOnTimestamps;
end


% wheelstart = findInSorted(spike2_data.analog_signal_time_vect, pupil_time(1));
% wheel_trace = spike2_data.wheelSpeed(wheelstart:500:end);
% wheel_trace = wheel_trace(1:max_pupil);
%corr(vid_energy, wheel_trace')

%% Grooming 

grooming = smooth(vid_energy, 20);          
%grooming_diff = diff(grooming);  


% figure();
% ax(1) = subplot(211);
% plot(pupil_time, vid_energy);
% ax(2) = subplot(212);
% plot(pupil_time, grooming);
% linkaxes(ax,'x');


thresh = 1.5; %was 0.3
[groomHighOn, groomHighOff] = cjb_changepoints(abs(grooming),pupil_time,thresh,1); 

% figure();
% plot(pupil_time, grooming); hold on;
% xline(groomHighOn, 'g')
% xline(groomHighOff, 'r');

[groomHighOn, groomHighOff] = mergeTimestamps(groomHighOn, groomHighOff, 1.2); %was 1

groomHighOn_temp=squeeze(groomHighOn(groomHighOff-groomHighOn>=1)); %was 1.5
groomHighOff_temp=squeeze(groomHighOff(groomHighOff-groomHighOn>=1)); %was 1.5

[groomHighOn_int2,groomHighOff_int2] = mergeTimestamps(groomHighOn_temp, groomHighOff_temp, 10);

if length(groomHighOn_int2) > 25
    disp(length(groomHighOn))
    thresh = 2; %was 2.5
    [groomHighOn, groomHighOff] = cjb_changepoints(abs(grooming),pupil_time,thresh,1);
    [groomHighOn, groomHighOff] = mergeTimestamps(groomHighOn, groomHighOff, 1.2); %was 1

    groomHighOn_temp=squeeze(groomHighOn(groomHighOff-groomHighOn>=1)); %was 1.5
    groomHighOff_temp=squeeze(groomHighOff(groomHighOff-groomHighOn>=1)); %was 1.5

    [groomHighOn_int2,groomHighOff_int2] = mergeTimestamps(groomHighOn_temp, groomHighOff_temp, 10);
end

groomHighOn_final=squeeze(groomHighOn_int2(groomHighOff_int2-groomHighOn_int2>=minGroomDur));
groomHighOff_final=squeeze(groomHighOff_int2(groomHighOff_int2-groomHighOn_int2>=minGroomDur));

% figure('Position', [50, 600, 1800, 200]);
% plot(pupil_time, vid_energy); hold on
% if ~isempty(groomHighOn_final)  
%     xline(groomHighOn_final, 'g')
%     xline(groomHighOff_final, 'r'); title(session, "Interpreter", 'none')
% end

if length(groomHighOn_final)>40
    warning('has more than 40 groom bouts, not including any')
    groomHighOn_int2 = [];
    groomHighOff_int2 = [];
    groomHighOn_final = [];
    groomHighOff_final = [];
end


%% og grooming
% grooming = smooth(vid_energy, 10);       
% %groom_thresh = quantile(grooming, 0.95);
% disp(['groom sum is : ' num2str(sum(grooming))])
% 
% % if sum(grooming) < 0
% %     thresh = 2;
% % else
% %     thresh = 3; %this may have been at 3 instead of 3.5
% % end
% 
% thresh = 2;
% 
% [groomHighOn, groomHighOff] = cjb_changepoints(grooming,spike2_data.pupilFrameOnTimestamps,thresh,1);
% 
% % [groomHighOn_int,groomHighOff_int] = mergeTimestamps(groomHighOn, groomHighOff, 2);
% % 
% % groomHighOn_int1=squeeze(groomHighOn_int(groomHighOff_int-groomHighOn_int>1)); %this had been at 1
% % groomHighOff_int1=squeeze(groomHighOff_int(groomHighOff_int-groomHighOn_int>1));
% 
% [groomHighOn_int2,groomHighOff_int2] = mergeTimestamps(groomHighOn, groomHighOff, 8);
% 
% % figure();
% % plot(pupil_time, grooming); hold on
% % xline(groomHighOn_int2, 'g')
% % xline(groomHighOff_int2, 'r')
% 
% groomHighOn_final=squeeze(groomHighOn_int2(groomHighOff_int2-groomHighOn_int2>=minGroomDur));
% groomHighOff_final=squeeze(groomHighOff_int2(groomHighOff_int2-groomHighOn_int2>=minGroomDur));
% 
% figure('Position', [100, 100, 900, 100]);
% plot(pupil_time, grooming); hold on
% if ~isempty(groomHighOn_final)
%     xline(groomHighOn_final, 'g')
%     xline(groomHighOff_final, 'r'); title(session, "Interpreter", 'none')
% end
% 
% if length(groomHighOn_final)>40
%     warning('has more than 40 groom bouts, not including any')
%     groomHighOn_int2 = [];
%     groomHighOff_int2 = [];
%     groomHighOn_final = [];
%     groomHighOff_final = [];
% end





%% Locomotion 


% remove running bouts before/after meso start/stop times
LocoOn = spike2_data.wheelOn >= imaging_time(1) & spike2_data.wheelOn <= imaging_time(end); % bouts before start
LocoOff = spike2_data.wheelOff >= imaging_time(1) & spike2_data.wheelOff <= imaging_time(end); % bouts after end 
wheelOn_int = spike2_data.wheelOn(LocoOn);
wheelOff_int = spike2_data.wheelOff(LocoOff);


% account for running when wheel recording started/ended
if wheelOn_int(1) > wheelOff_int(1) % account for animal running when wheel recording started
    wheelOff_int = wheelOff_int(2:end);
end

if wheelOff_int(end) < wheelOn_int(end) % account for animal running when wheel recording started
    wheelOff_int = [wheelOff_int; imaging_time(end)];
end 


% make sure no grooming is co-occuring
[notGroomOn, notGroomOff] = timestampsNot(groomHighOn_int2, groomHighOff_int2, 'StartTime', 0, 'EndTime', wheelTime(end));
[wheelOn_int2, wheelOff_int2] = timestampsAnd(wheelOn_int, wheelOff_int, notGroomOn, notGroomOff);



% figure();
% plot(spike2_data.analog_signal_time_vect, spike2_data.wheelSpeed); hold on
% xline(wheelOn_int2, 'g')
% xline(wheelOff_int2, 'r')



locoBoutTime_int = wheelOff_int2 - wheelOn_int2 > 0.25;
wheelOn_int2 = wheelOn_int2(locoBoutTime_int);
wheelOff_int2 = wheelOff_int2(locoBoutTime_int);

% merge loco bouts within 3s of eachother
[wheelOn_final, wheelOff_final] = mergeTimestamps(wheelOn_int2, wheelOff_int2, 3);


% figure();
% plot(spike2_data.analog_signal_time_vect, spike2_data.wheelSpeed); hold on
% xline(wheelOn_final, 'g')
% xline(wheelOff_final, 'r')

locoBoutTime = wheelOff_final - wheelOn_final > minRunDur;
wheelOn_final = wheelOn_final(locoBoutTime);
wheelOff_final = wheelOff_final(locoBoutTime);


%% Sit (Quiescience)


sitOn = spike2_data.wheelOff;
sitOff = spike2_data.wheelOn;

SitOn = sitOn >= imaging_time(1) & sitOn <= imaging_time(end);
SitOff = sitOff >= imaging_time(1) & sitOff <= imaging_time(end);
sitOn_int = sitOn(SitOn);
sitOff_int = sitOff(SitOff);


% account for sitting when wheel recording started/ended
if sitOn_int(1) > sitOff_int(1) % account for animal running when wheel recording started
    sitOn_int = [imaging_time(1); sitOn_int];
end
 
if sitOff_int(end) < sitOn_int(end) % account for animal running when wheel recording started
    sitOff_int = [sitOff_int; imaging_time(end)];
end


[notGroomOn, notGroomOff] = timestampsNot(groomHighOn_int2, groomHighOff_int2, 'StartTime', 0, 'EndTime', wheelTime(end));
[sitOn_final, sitOff_final] = timestampsAnd(sitOn_int, sitOff_int, notGroomOn, notGroomOff);

% need to threshold for total bout time after removing timepoints that overlap with grooming
sitBoutTime = sitOff_final - sitOn_final > minSitDur;
sitOn_final = sitOn_final(sitBoutTime);
sitOff_final = sitOff_final(sitBoutTime);


if isempty(sitOn_final)
    warning('no sitting bouts detected')
    states.locoOn = wheelOn_final' + TimeSinceLocOn;
    states.locoOff = wheelOff_final' - TimeBeforeLocOff;
    states.sitOn = [];
    states.sitOff = [];
    states.faceHighSitOn = [];
    states.faceHighSitOff = [];
    states.faceLowSitOn = [];
    states.faceLowSitOff = [];
    states.groomHighOn = groomHighOn_int2;
    states.groomHighOff = groomHighOff_int2; 
    return
else
    sitOn_final = sitOn_final + TimeSinceSitOn;
    sitOff_final = sitOff_final - TimeBeforeSitOff; 
end



%% FACE MAP



% gets time points within sitting bouts 
pupilTime_Idx=cell(1,length(sitOn_final));
for st=1:length(sitOn_final)
    pupilTime_Idx{st}=find(pupil_time > sitOn_final(st) & pupil_time < sitOff_final(st));
end
pupilTime_quiescence=cell2mat(pupilTime_Idx');

% normalized face frameMotion during quiescent bouts 
face_quiescence = vid_energy(pupilTime_quiescence);

zThresh_high=quantile(face_quiescence,0.60);
zThresh_low=quantile(face_quiescence,0.40);


% High Face 
[faceHighOn,faceHighOff] = cjb_changepoints(vid_energy, pupil_time, zThresh_high, 1); 

[faceHighSitOn,faceHighSitOff] = timestampsAnd(faceHighOn,faceHighOff,sitOn_final,sitOff_final);
faceHighOn_inter=squeeze(faceHighSitOn(faceHighSitOff-faceHighSitOn>minStateDur));
faceHighOff_inter=squeeze(faceHighSitOff(faceHighSitOff-faceHighSitOn>minStateDur));

% figure();
% plot(pupil_time, motion_energy); 
% yline(zThresh_high)
% xline(faceHighOn_inter, 'g')
% xline(faceHighOff_inter, 'r')

% [alignedData, ~] = averageMovieToEvent(vid_energy', pupil_time, faceHighOn_inter, 5, 15);
% mean_whiskOn = squeeze(mean(alignedData, 1));
% figure();plot(-4:15, mean_whiskOn); hold on
% xline(0); title('High Face Onset'); xlabel('time (s)');ylabel ('zscore(facial motion)')



% Low Face
[faceLowOn,faceLowOff] = cjb_changepoints(-vid_energy, pupil_time, -zThresh_low, 1); 

[faceLowSitOn,faceLowSitOff] = timestampsAnd(faceLowOn,faceLowOff,sitOn_final,sitOff_final);
faceLowSitOn_inter = squeeze(faceLowSitOn(faceLowSitOff-faceLowSitOn>minStateDur));
faceLowSitOff_inter = squeeze(faceLowSitOff(faceLowSitOff-faceLowSitOn>minStateDur));
 

% figure();
% plot(pupil_time, motion_energy); 
% yline(zThresh_low)
% xline(faceLowSitOn_inter, 'g')
% xline(faceLowSitOff_inter, 'r')
% 
% 
% [alignedData, ~] = averageMovieToEvent(vid_energy', pupil_time, faceLowSitOn_inter, 5, 15);
% mean_quiOn = squeeze(mean(alignedData, 1));
% figure();plot(-4:15, mean_quiOn); hold on
% xline(0);title('Low Face Onset'); xlabel('time (s)'); ylabel ('zscore(facial motion)')
% 


%% add them all to states struct 

% add to states data
states.locoOn = wheelOn_final' + TimeSinceLocOn;
states.locoOff = wheelOff_final' - TimeBeforeLocOff;

states.sitOn = sitOn_final; 
states.sitOff = sitOff_final;

states.faceHighSitOn = faceHighOn_inter ; 
states.faceHighSitOff = faceHighOff_inter ;

states.faceLowSitOn = faceLowSitOn_inter ;
states.faceLowSitOff = faceLowSitOff_inter ;

% states.faceLowLowSitOn = faceLowLowSitOn_inter ;
% states.faceLowLowSitOff = faceLowLowSitOff_inter ;

states.groomHighOn = groomHighOn_final;
states.groomHighOff = groomHighOff_final;


end



