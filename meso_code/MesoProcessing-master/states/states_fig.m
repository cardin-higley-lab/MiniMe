function fig = states_fig(states, spike2_data, vid_energy, animal, session)


%%


% if isfield(spike2_data, 'blueOnTimestamps')
%     imaging_time = spike2_data.blueOnTimestamps;
% else
%     imaging_time = spike2_data.greenOnTimestamps;
% end
% 
% imaging_start_time = imaging_time(1);
% 
% 
% behavior_start_ind = findClosestDouble(spike2_data.pupilFrameOnTimestamps, imaging_start_time);
% pupil_time = spike2_data.pupilFrameOnTimestamps(behavior_start_ind:end);
% vid_energy = vid_energy(behavior_start_ind:end);
% 
% maxL = min(length(pupil_time), length(vid_energy));
% pupil_time = pupil_time(1:maxL);
% vid_energy = vid_energy(1:maxL);


%wheelTime = (1:length(spike2_data.wheelSpeed))/fsspike2;

pupil_time = spike2_data.pupilFrameOnTimestamps;
maxL = min(length(pupil_time), length(vid_energy));
pupil_time = pupil_time(1:maxL);
vid_energy = vid_energy(1:maxL);

%%


% Visualize the States

if ~isempty(vid_energy)
    fig = figure();
    
    if ~isempty(spike2_data.wheelSpeed)
        ax(1) = subplot(511);
        plot(spike2_data.analog_signal_time_vect, spike2_data.wheelSpeed);
        ylabel('speed (cm/s)');
        if ~isempty(states.locoOn)
            xline(states.locoOn, 'g');
            xline(states.locoOff, 'r');
        end
        title('Locomotions bouts')
        
        ax(2) = subplot(512);
        plot(spike2_data.analog_signal_time_vect, spike2_data.wheelSpeed);
        ylabel('speed (cm/s)');
        if ~isempty(states.sitOn)
            xline(states.sitOn, 'g');
            xline(states.sitOff, 'r');
        end
        title('Sitting bouts')
    

    end
    
    ax(3) = subplot(513);
    plot(pupil_time, vid_energy);
    if ~isempty(states.groomHighOn)
        xline(states.groomHighOn, 'g')
        xline(states.groomHighOff, 'r')
    end
    ylabel('grooming vid_energy1')
    title('Grooming bouts')
    
    ax(4) = subplot(514);
    plot(pupil_time, vid_energy);
    ylabel('face vid_energy1');
    if ~isempty(states.faceHighSitOn)
        xline(states.faceHighSitOn, 'g');
        xline(states.faceHighSitOff, 'r');
    end
    %hold on; plot(spike2_data.pupilFrameOnTimestamps(1:maxlength), boolResult(1:maxlength)*1000)
    title('Face high & sitting')
    
    ax(5) = subplot(515);
    plot(pupil_time, vid_energy);
    if ~isempty(states.faceLowSitOn)
        xline(states.faceLowSitOn, 'g')
        xline(states.faceLowSitOff, 'r')
    end
    ylabel('face vid_energy1')
    title('Face Low & sitting')
    linkaxes(ax,'x');
    
    sgtitle(['session: ' , animal ', ', session], 'Interpreter', 'none');

elseif isempty(vid_energy)
    fig = figure();
    
    if ~isempty(spike2_data.wheelSpeed)
        ax(1) = subplot(211);
        plot(spike2_data.analog_signal_time_vect, spike2_data.wheelSpeed);
        ylabel('speed (cm/s)');
        if ~isempty(states.locoOn)
            xline(states.locoOn, 'g');
            xline(states.locoOff, 'r');
        end
        title('Locomotions bouts')
        
        ax(2) = subplot(212);
        plot(spike2_data.analog_signal_time_vect, spike2_data.wheelSpeed);
        ylabel('speed (cm/s)');
        if ~isempty(states.sitOn)
            xline(states.sitOn, 'g');
            xline(states.sitOff, 'r');
        end
        title('Sitting bouts')
    
  
    end
    
    linkaxes(ax,'x');
    sgtitle(['session: ' , animal ', ', session], 'Interpreter', 'none');
end