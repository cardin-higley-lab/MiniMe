control_Ts = [controlDeviant1TimestampsOn; controlDeviant2TimestampsOn];
control_amp_face = nan(length(control_Ts), 2);
for i = 1:length(control_Ts)

    ind = findInSorted(mesoTimestamps,control_Ts(i));
    for groombout = 1:length(states.groomHighOn)
        if  mesoTimestamps(ind) > states.groomHighOn(groombout) &&  mesoTimestamps(ind) < states.groomHighOff(groombout)
            disp('true');
            continue;
        end
    end
    v1_amplitude_before = sum(v1_data(ind-9:ind));
    v1_amplitude_after = sum(v1_data(ind:ind+9));
    v1_amplitude = v1_amplitude_after - v1_amplitude_before;
    face_val = mean(face_Norm(ind-4:ind));
    control_amp_face(i, :) = [v1_amplitude, face_val];
    
end

redundant_Ts = [redundant2TimestampsOn;];% redundant1TimestampsOn; 
redundant_amp_face = nan(length(redundant_Ts), 2);
for i = 1:length(redundant_Ts)

    ind = findInSorted(mesoTimestamps,redundant_Ts(i));
    for groombout = 1:length(states.groomHighOn)
        if mesoTimestamps(ind) > states.groomHighOn(groombout) && mesoTimestamps(ind) < states.groomHighOff(groombout)
            disp('true');
            continue;
        end
    end
    v1_amplitude_before = sum(v1_data(ind-9:ind));
    v1_amplitude_after = sum(v1_data(ind:ind+9));
    v1_amplitude = v1_amplitude_after - v1_amplitude_before;
    face_val = mean(face_Norm(ind-4:ind));
    redundant_amp_face(i, :) = [v1_amplitude, face_val];
    
end

deviant_Ts = [deviant1TimestampsOn; ];%deviant2TimestampsOn 
deviant_amp_face = nan(length(deviant_Ts), 2);
for i = 1:length(deviant_Ts)

    ind = findInSorted(mesoTimestamps,deviant_Ts(i));
    for groombout = 1:length(states.groomHighOn)
        if  mesoTimestamps(ind) > states.groomHighOn(groombout) &&  mesoTimestamps(ind) < states.groomHighOff(groombout)
            disp('true');
            continue;
        end
    end
    v1_amplitude_before = sum(v1_data(ind-9:ind));
    v1_amplitude_after = sum(v1_data(ind:ind+9));
    v1_amplitude = v1_amplitude_after - v1_amplitude_before;
    face_val = mean(face_Norm(ind-4:ind));
    deviant_amp_face(i, :) = [v1_amplitude, face_val];
    
end



%bar plot 
x=categorical({'Control', 'Redundant', 'Deviant'});
x=reordercats(x, {'Control', 'Redundant', 'Deviant'});
y_means = [mean(control_amp_face(:,1)); mean(redundant_amp_face(:,1)); mean(deviant_amp_face(:,1))];
barplot = figure; hold on; box off
bar(x,y_means); 
ylabel('AUC 1s following vis stim') %,'FontSize',18
saveas(barplot, 'CvsRvsD.png')
