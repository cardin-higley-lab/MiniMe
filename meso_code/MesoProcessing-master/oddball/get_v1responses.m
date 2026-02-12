function output = get_v1responses(timestamps, mesoTimestamps, states, v1_data, face_Norm)

output = nan(length(timestamps), 3);
for i = 1:length(timestamps)

    ind = findInSorted(mesoTimestamps,timestamps(i));
    for groombout = 1:length(states.groomHighOn)
        if mesoTimestamps(ind) > states.groomHighOn(groombout) && mesoTimestamps(ind) < states.groomHighOff(groombout)
            disp('true');
            continue;
        end
    end
    v1_amplitude_before = mean(v1_data(ind-20:ind));
    v1_amplitude_after = mean(v1_data(ind:ind+9));
    v1_amplitude = v1_amplitude_after - v1_amplitude_before;
    face_val = mean(face_Norm(ind-4:ind));
    output(i, :) = [v1_amplitude, face_val, i];
    
end