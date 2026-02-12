function processedSpike2Data = processSpike2Data(spike2Data, vis_isi)
%%processSpike2 gets relevant info from spike2 files 
% Data is extracted from spike2 files via a dedicated windows specific
% library, then processed according to channel type. The majority of analog
% signals are decomposed into states via kmeans clustering with the
% analogSigToTimestamps function. Some signals are unaltered, such as
% EEG. The wheel signal is transformed into wheel speed, and periods of
% running are identified with the extract_wheel_signal function.

spike2Fields = fields(spike2Data);
processedSpike2Data = struct;
for i = 1:length(spike2Fields)
        currDataBundle = spike2Data.(spike2Fields{i});
        timestamps = cumsum(ones(size(currDataBundle.data))/currDataBundle.samplingRate)+double(currDataBundle.acquisitionOffset);
        if contains(spike2Fields{i},'GRLED') || contains(spike2Fields{i},'Green') || contains(spike2Fields{i},'GreenLED') || contains(spike2Fields{i},'Green LED') 
            [protoGreenOnTimestamps,protoGreenOffTimestamps]=analogSigToTimestamps(spike2Data.(spike2Fields{i}).data,timestamps);
            processedSpike2Data.greenOnTimestamps = protoGreenOnTimestamps; % truncate to match dFoF (calculated elsewhere)
            processedSpike2Data.greenOffTimestamps = protoGreenOffTimestamps;
%             figure(); plot(spike2Data.(spike2Fields{i}).data); title("green")
%             processedSpike2Data.temp1 = spike2Data.(spike2Fields{i}).data;

        elseif contains(spike2Fields{i},'BLLED') || contains(spike2Fields{i},'BLUELED')
            [protoBlueOnTimestamps,protoBlueOffTimestamps]=analogSigToTimestamps(spike2Data.(spike2Fields{i}).data,timestamps);
            processedSpike2Data.blueOnTimestamps = protoBlueOnTimestamps; % truncate to match dFoF (calculated elsewhere)
            processedSpike2Data.blueOffTimestamps = protoBlueOffTimestamps;

        elseif  contains(spike2Fields{i},'BACKSCAT') || contains(spike2Fields{i},'bs') || contains(spike2Fields{i},'BSLED') || contains(spike2Fields{i},'BS') 
            [protoGreenOnTimestamps,protoGreenOffTimestamps]=analogSigToTimestamps(spike2Data.(spike2Fields{i}).data,timestamps);
            processedSpike2Data.greenOnTimestamps = protoGreenOnTimestamps; % truncate to match dFoF (calculated elsewhere)
            processedSpike2Data.greenOffTimestamps = protoGreenOffTimestamps;
            %figure(); plot(spike2Data.(spike2Fields{i}).data); title("backscat")
            %disp(spike2Data.(spike2Fields{i}).data(1:10, 1))
%             processedSpike2Data.temp2 = spike2Data.(spike2Fields{i}).data;

        elseif contains(spike2Fields{i},'UV') || contains(spike2Fields{i},'UVLED')
            [protoUvOnTimestamps,protoUvOffTimestamps]=analogSigToTimestamps(spike2Data.(spike2Fields{i}).data,timestamps);
             processedSpike2Data.uvOnTimestamps = protoUvOnTimestamps;
             processedSpike2Data.uvOffTimestamps = protoUvOffTimestamps;
        elseif contains(spike2Fields{i},'MesoCam')
            [protoMesoFrameOnTimestamps,protoMesoFrameOffTimestamps]=analogSigToTimestamps(spike2Data.(spike2Fields{i}).data,timestamps);
            processedSpike2Data.mesoFrameOnTimestamps = protoMesoFrameOnTimestamps;
            processedSpike2Data.mesoFrameOffTimestamps = protoMesoFrameOffTimestamps;
        elseif contains(spike2Fields{i},'Vis')
            %[processedSpike2Data.diodeOnTimestamps,processedSpike2Data.diodeOffTimestamps] = analogSigToTimestamps(spike2Data.(spike2Fields{i}).data,timestamps);

            diode_data = spike2Data.(spike2Fields{i}).data;
            if min(diode_data) < -0.265
                above_thresh = find(diode_data > -.26);
                if ~isempty(above_thresh)
                    diode_data(above_thresh) = mean(diode_data, 'omitnan');
                end
                [onDiode,offDiode] = analogSigToTimestamps(diode_data,timestamps);
                [mergeOn,mergeOff] = mergeTimestamps(onDiode,offDiode,vis_isi);
                processedSpike2Data.diodeOnTimestamps = mergeOn;
                processedSpike2Data.diodeOffTimestamps = mergeOff;
                time = (0:(length(spike2Data.(spike2Fields{i}).data)-1))/5000;
                visOn = zeros(1, length(time));
                for ii = 1:length(mergeOn)
                    onT = findInSorted(time, mergeOn(ii));
                    offT = findInSorted(time, mergeOff(ii));
                    visOn(onT:offT) = 1;
                end
                %processedSpike2Data.vis_trace = onDiode; warning('mergeTimestamps off')
                processedSpike2Data.vis_trace = visOn;
            end 


        elseif contains(spike2Fields{i},'wheel')
            [processedSpike2Data.wheelSpeed,wheel_on,wheel_off] = cjb_extract_wheel_signal(spike2Data.(spike2Fields{i}).data);
            processedSpike2Data.wheelOn = timestamps(wheel_on);
            processedSpike2Data.wheelOff = timestamps(wheel_off);
            processedSpike2Data.analog_signal_time_vect = (0:(length(processedSpike2Data.wheelSpeed)-1))/5000;

        elseif contains(spike2Fields{i},'piezo')
            processedSpike2Data.piezo = spike2Data.(spike2Fields{i}).data;
            processedSpike2Data.wheelSpeed = [];
            processedSpike2Data.wheelOn = [];
            processedSpike2Data.wheelOff = [];
            processedSpike2Data.analog_signal_time_vect = (0:(length(spike2Data.(spike2Fields{i}).data)-1))/5000;

        elseif contains(spike2Fields{i},'Air') || contains(spike2Fields{i},'airpuff')
            [processedSpike2Data.airPuffOnTimestamps,processedSpike2Data.airPuffOffTimestamps]=analogSigToTimestamps(spike2Data.(spike2Fields{i}).data,timestamps);
        elseif contains(spike2Fields{i},'pupil')
            [processedSpike2Data.pupilFrameOnTimestamps,processedSpike2Data.pupilFrameOffTimestamps]=analogSigToTimestamps(spike2Data.(spike2Fields{i}).data,timestamps);
        elseif contains(spike2Fields{i},'Water') || contains(spike2Fields{i},'water') 
            [processedSpike2Data.waterOnTimestamps,processedSpike2Data.waterOffTimestamps]=analogSigToTimestamps(spike2Data.(spike2Fields{i}).data,timestamps);
        elseif contains(spike2Fields{i},'body') || contains(spike2Fields{i},'Body') || contains(spike2Fields{i},'BODY')
            [processedSpike2Data.bodyFrameOnTimestamps,processedSpike2Data.bodyFrameOffTimestamps]=analogSigToTimestamps(spike2Data.(spike2Fields{i}).data,timestamps);
        elseif contains(spike2Fields{i},'LFP') || contains(spike2Fields{i},'ECoG')
            processedSpike2Data.ECoG = spike2Data.(spike2Fields{i}).data;
        else
            [processedSpike2Data.([spike2Fields{i} 'OnTimestamps']), processedSpike2Data.([spike2Fields{i} 'OffTimestamps'])]=analogSigToTimestamps(spike2Data.(spike2Fields{i}).data,timestamps);
            warning(['Unindetified channel name: ' spike2Fields{i} ', treating as timestamp data']);
        end
end
end
