function states = locoTimestamps(spike2_data, state_duration)
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

minStateDur = state_duration; %making all states need to be at least 5 sec to count 
TimeSinceLocOn=3;%for locomotion state, minimum time since locomotion onset
TimeBeforeLocOff=3;%for locomotion state, minimum time before locomotion offset
TimeSinceSitOn=10;%for quiescence state, minimum time since quiescence onset
TimeBeforeSitOff=10;%for quiescence state, minimum time before quiescence offset
fsspike2=5000;

minSitDur=minStateDur+TimeSinceSitOn+TimeBeforeSitOff; 
minRunDur=minStateDur+TimeSinceLocOn+TimeBeforeLocOff;
imaging_time = spike2_data.blueOnTimestamps;
wheelTime = (1:length(spike2_data.wheelSpeed))/fsspike2;

%% Locomotion 

if isfield(spike2_data, 'piezo')
    disp('piezo used, no locomotion')
    wheelOn_final = [];
    wheelOff_final = [];
elseif  max(spike2_data.wheelSpeed) < 1
    wheelOn_final = [];
    wheelOff_final = [];
else

    wheel_on = spike2_data.wheelOn; 
    wheel_off = spike2_data.wheelOff; 

    % remove running bouts before/after meso start/stop times 
    firstLocoOn = wheel_on >= imaging_time(1); % first wheel On should not be before start of imaging 
    lastLocoOff = wheel_off <= imaging_time(end); % last wheel Off should not be after end of imaging 
    wheelOn_int = wheel_on(firstLocoOn);
    wheelOff_int = wheel_off(lastLocoOff);

    % account for running when wheel recording started/ended
    if wheelOn_int(1) > wheelOff_int(1) % account for animal running when wheel recording started
        wheelOn_int = [imaging_time(1); wheelOn_int];
    end
    
    if wheelOff_int(end) < wheelOn_int(end) % account for animal running when wheel recording started
        wheelOff_int = [wheelOff_int; imaging_time(end)];
    end
   
    % if animal was running when started/ended, crop bout by 1 sec
    for whe=1:length(wheelOn_int)
        if wheelOn_int(whe)<imaging_time(end) && wheelOff_int(whe)>imaging_time(end)
            wheelOff_int(whe)=imaging_time(end)-1;%if locomotion starts before end of imaging but continues after, only extract state until imaging time end minus a second
        end
        
        if wheelOn_int(whe)<imaging_time(1) && wheelOff_int(whe)>imaging_time(1)
            wheelOn_int(whe)=imaging_time(1)+1;%if locomotion starts before start of imaging but continues after, only extract state from imaging time start plus a second
        end
    end
    
    wheelOn_int=(wheelOn_int(:))'; 
    wheelOff_int=(wheelOff_int(:))';
    
    
    % remove running bouts where speed does not go above 3 cm/s
    wheelOn_int2 = zeros(1, length(wheelOn_int));
    wheelOff_int2 = zeros(1, length(wheelOff_int));
    max_speeds = zeros(1, length(wheelOff_int));
    
    for loco_bout = 1:length(wheelOn_int)
        ind1 = findInSorted(spike2_data.analog_signal_time_vect, wheelOn_int(loco_bout));
        ind2 = findInSorted(spike2_data.analog_signal_time_vect, wheelOff_int(loco_bout));
        speed_max = max(spike2_data.wheelSpeed(ind1:ind2));
        max_speeds(loco_bout) = speed_max;
        if speed_max > 3
            wheelOn_int2(loco_bout) = wheelOn_int(loco_bout); 
            wheelOff_int2(loco_bout) = wheelOff_int(loco_bout);
        else
            wheelOn_int2(loco_bout) = NaN;
            wheelOff_int2(loco_bout) = NaN;
        end
    end
    wheelOn_int3 = wheelOn_int2(~isnan(wheelOn_int2));
    wheelOff_int3 = wheelOff_int2(~isnan(wheelOff_int2));

    % merge loco bouts within 2s of eachother 
    wheelOn_int4 = wheelOn_int3;
    wheelOff_int4 = wheelOff_int3;
    for i = 1:(length(wheelOn_int3)-1)
        if wheelOn_int3(i+1) < (wheelOff_int3(i) + 2)
            wheelOn_int4(i+1) = nan;
            wheelOff_int4(i) = nan;
        end
    end
    wheelOn_int4 = wheelOn_int4(~isnan(wheelOn_int4));
    wheelOff_int4 = wheelOff_int4(~isnan(wheelOff_int4));

    locoBoutTime = wheelOff_int4 - wheelOn_int4 > minRunDur; 
    wheelOn_final = wheelOn_int4(locoBoutTime);
    wheelOff_final = wheelOff_int4(locoBoutTime);

end


%% Sit (Quiescience)

if isfield(spike2_data, 'piezo')
    sitOn_final = imaging_time(1);
    sitOff_final = imaging_time(end);
elseif max(spike2_data.wheelSpeed) < 1
    sitOn_final = imaging_time(1);
    sitOff_final = imaging_time(end);
else
    sitOn = spike2_data.wheelOff;
    sitOff = spike2_data.wheelOn; 
    
    firstSitOn = sitOn >= imaging_time(1);
    lastSitOff = sitOff <= imaging_time(end);
    sitOn_int = sitOn(firstSitOn);
    sitOff_int = sitOff(lastSitOff);

    % account for sitting when wheel recording started/ended
    if sitOn_int(1) > sitOff_int(1) % account for animal running when wheel recording started
        sitOn_int = [imaging_time(1); sitOn_int];
    end
    
    if sitOff_int(end) < sitOn_int(end) % account for animal running when wheel recording started
        sitOff_int = [sitOff_int; imaging_time(end)];
    end
    
    % need to threshold for total bout time after removing timepoints that overlap with grooming 
    sitBoutTime = sitOff_int - sitOn_int > minSitDur; 
    sitOn_final = sitOn_int(sitBoutTime);
    sitOff_final = sitOff_int(sitBoutTime);

end


%% add them all to states struct 

% add to states data
states.locoOn = wheelOn_final' + TimeSinceLocOn;
states.locoOff = wheelOff_final' - TimeBeforeLocOff;

states.sitOn = sitOn_final + TimeSinceSitOn;
states.sitOff = sitOff_final - TimeBeforeSitOff; 


end



