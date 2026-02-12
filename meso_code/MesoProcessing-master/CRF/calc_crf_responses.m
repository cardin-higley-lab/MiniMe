function [vis_responses, loco] = calc_crf_responses(contrasts, parcels_time_trace, spike2_data, states)

% calculates crf responses sitting vs running for an individual animal based on contrast order and behavior

[numParcels, ~] = size(parcels_time_trace);

if isfield(spike2_data, 'blueOnTimestamps')
    mesoOnTimestamps = spike2_data.blueOnTimestamps;
else
    mesoOnTimestamps = spike2_data.greenOnTimestamps;
end
wheelTime = (1:length(spike2_data.wheelSpeed))/5000;


%% clean up wheel data

wheel_on = spike2_data.wheelOn;
wheel_off = spike2_data.wheelOff;

% remove running bouts before/after meso start/stop times
LocoOn = wheel_on >= mesoOnTimestamps(1) & wheel_on <= mesoOnTimestamps(end); % first wheel On should not be before start of imaging
LocoOff = wheel_off >= mesoOnTimestamps(1) & wheel_off <= mesoOnTimestamps(end); % last wheel Off should not be after end of imaging
wheelOn_int = wheel_on(LocoOn);
wheelOff_int = wheel_off(LocoOff);

% account for running when wheel recording started/ended
if wheelOn_int(1) > wheelOff_int(1) % account for animal running when wheel recording started
    wheelOn_int = [mesoOnTimestamps(1); wheelOn_int];
end

if wheelOff_int(end) < wheelOn_int(end) % account for animal running when wheel recording started
    wheelOff_int = [wheelOff_int; mesoOnTimestamps(end)];
end

% if animal was running when started/ended, crop bout by 1 sec
for whe=1:length(wheelOn_int)
    if wheelOn_int(whe)<mesoOnTimestamps(end) && wheelOff_int(whe)>mesoOnTimestamps(end)
        wheelOff_int(whe)=mesoOnTimestamps(end)-1;%if locomotion starts before end of imaging but continues after, only extract state until imaging time end minus a second
    end

    if wheelOn_int(whe)<mesoOnTimestamps(1) && wheelOff_int(whe)>mesoOnTimestamps(1)
        wheelOn_int(whe)=mesoOnTimestamps(1)+1;%if locomotion starts before start of imaging but continues after, only extract state from imaging time start plus a second
    end
end


wheelOn_int=(wheelOn_int(:))';
wheelOff_int=(wheelOff_int(:))';


[notGroomOn, notGroomOff] = timestampsNot(states.groomHighOn, states.groomHighOn, 'StartTime', 0, 'EndTime', wheelTime(end));
[wheelOn_final, wheelOff_final] = timestampsAnd(wheelOn_int, wheelOff_int, notGroomOn, notGroomOff);


% find stims that count as locomotion (includes if running ended 0.5s before stim onset or started 1s after stim offset 
stim_in_loco = false(length(mesoOnTimestamps),1);
for ii = 1:length(wheelOn_final)
    stim_in_loco(mesoOnTimestamps >= wheelOn_final(ii)-0.5 & mesoOnTimestamps <= wheelOff_final(ii) + 1)  = true;
end




%% get crf stuff now...

contrast_vals = unique(contrasts);
numContrasts = length(contrast_vals);
numPres = length(contrasts)/numContrasts;
contrast_inds = nan(numContrasts, numPres);
contrast_times = nan(numContrasts, numPres);
vis_responses = nan(numParcels, numContrasts, numPres);
loco = nan(numContrasts, numPres);

for parcels = 1:numParcels
    
    data = mean(parcels_time_trace(parcels,:), 1, 'omitnan');

    for i = 1:numContrasts
        contrast = contrast_vals(i);
        inds = find(contrasts == contrast);
        contrast_inds(i, :) = inds;
        times = spike2_data.diodeOnTimestamps(inds);
        contrast_times(i, :) = times;
    
        vis_response = nan(numPres, 1);
        loco_val = nan(numPres, 1);
        
    
        for ii = 1:numPres
            stim_time = times(ii);
            
            % first looking for individual vis responses 
            ind = findInSorted(mesoOnTimestamps,stim_time);
            v1_before = mean(data(ind-9:ind));
            v1_after = mean(data(ind:ind+20));
            vis_response(ii, 1) = v1_after/v1_before;
            
            % also need to see if loco or not loco 
            loco_val(ii, 1) = stim_in_loco(ind);
               
        end
        
        loco(i, :) = loco_val;
        parcel_responses(i, :) = vis_response;
    end
    vis_responses(parcels, :, :) = parcel_responses;

end


 

end


