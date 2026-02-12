function channel_data = get_relevant_spike2_data(dataSmrxFile,params)
%%get_relevant_spike2_data gets relevant info from spike2 files 
% Data is extracted from spike2 files via a dedicated windows specific
% library, then processed according to channel type. The majority of analog
% signals are decomposed into states via kmeans clustering with the
% analog_sig_to_timestamps function. Some signals are unaltered, such as
% EEG. The wheel signal is transformed into wheel speed, and periods of
% running are identified with the extract_wheel_signal function.
 

CEDS64LoadLib(params.cedpath) % load in the library to read spike2

disp(['loading ' dataSmrxFile]);

 if ~exist(dataSmrxFile,'file') % check if file exists
        error('Could not find file');
 else 
    fhand = CEDS64Open(dataSmrxFile); % open if it exists
    if fhand == -1 % check if it opened correctly
        error('Could not open file');
    end
 end

 
iChanNum = min(CEDS64MaxChan(fhand),length(fields(params.channels))); % get number of channels

dsec=CEDS64TimeBase(fhand); % seconds per sample

maxTimeTicks = CEDS64ChanMaxTime(fhand,1);
data=nan(maxTimeTicks./20+1,iChanNum);

% get analog data from each channel
readF = zeros([iChanNum,1]).*nan;
for ichan=1:iChanNum
    %file name, channel num, max num of points, start and end time in ticks
    [fRead,fVals,fTime] = CEDS64ReadWaveF(fhand,ichan,maxTimeTicks,0,maxTimeTicks);
    if fRead > 0
        data(fTime+1:fRead+fTime,ichan)=fVals; % save read data
        readF(ichan) = fRead; % store number of samples read for each channel
    end
end

%%
if length(unique(readF))~=1
    warning('Not all channels have the same length!');
end

data=data(1:readF(1)+fTime,:); % not sure what this does
timestamp=(1:20:maxTimeTicks)'*dsec; % create timestamps for each sample (not sure why downsampled by 20, maybe clock to sampling freq factor?

totallength=min(size(data,1),length(timestamp));
data=data(1:totallength,:); % trim some extra samples
channel_data.analog_signal_time_vect = timestamp(1:totallength);


names = fieldnames(params.channels); % get names of each channel to iterate through
% Here, we use the analog data in a channel specific manner. Most signals
% are composed of 2 states, thus we employ kmeans clustering through the
% analog_to_timestamps function to convert to timestamps, which we then
% save. Some functions like the wheel may have a more specific function,
% see below, and others are not altered (like EEG). The blue, UV, and
% frameticks/red meso channels are also truncated to match the truncation
% from filtering, allowing them to match up with dFoF (should really be
% end-params.deterend.filtLen/2-1, but dFoF is not truncated that way so it
% is not here)

for ni = 1:length(names)
    switch names{ni}
        case 'MOVING_PHOTO_DIODE'
            channel_data.movingdiode=data(:,params.channels.MOVING_PHOTO_DIODE);
        case 'BLUE'
            [protoBlueOnTimestamps,protoBlueOffTimestamps]=analog_sig_to_timestamps(data(:,params.channels.BLUE),channel_data.analog_signal_time_vect);
            channel_data.blueOnTimestamps = protoBlueOnTimestamps(params.deterend.filtLen/2+1:end-params.deterend.filtLen/2); % truncate to match dFoF (calculated elsewhere)
            channel_data.blueOffTimestamps = protoBlueOffTimestamps(params.deterend.filtLen/2+1:end-params.deterend.filtLen/2);
        case 'UV'
            [protoUvOnTimestamps,protoUvOffTimestamps]=analog_sig_to_timestamps(data(:,params.channels.UV),channel_data.analog_signal_time_vect);
             channel_data.uvOnTimestamps = protoUvOnTimestamps(params.deterend.filtLen/2+1:end-params.deterend.filtLen/2);
             channel_data.uvOffTimestamps = protoUvOffTimestamps(params.deterend.filtLen/2+1:end-params.deterend.filtLen/2);
        case {'FRAMETICKS' 'RED_MESO'}
            [protoMesoFrameOnTimestamps,protoMesoFrameOffTimestamps]=analog_sig_to_timestamps(data(:,params.channels.FRAMETICKS),channel_data.analog_signal_time_vect);
            channel_data.mesoFrameOnTimestamps = protoMesoFrameOnTimestamps(params.deterend.filtLen+1:end-params.deterend.filtLen);
            channel_data.mesoFrameOffTimestamps = protoMesoFrameOffTimestamps(params.deterend.filtLen+1:end-params.deterend.filtLen);
        case 'PHOTO_DIODE'
            [channel_data.diodeOnTimestamps,channel_data.diodeOffTimestamps]=analog_sig_to_timestamps(data(:,params.channels.PHOTO_DIODE),channel_data.analog_signal_time_vect);
        case 'WHEEL'
            [channel_data.wheelSpeed,protoWheelOn,protoWheelOff] = cjb_extract_wheel_signal(data(:, params.channels.WHEEL));
            channel_data.wheelOn = channel_data.analog_signal_time_vect(protoWheelOn); % wheel on and off in indices, so convert to seconds
            channel_data.wheelOff = channel_data.analog_signal_time_vect(protoWheelOff);
        case 'AIR_PUFF'
            [channel_data.airPuffOnTimestamps,channel_data.airPuffOffTimestamps]=analog_sig_to_timestamps(data(:,params.channels.AIR_PUFF),channel_data.analog_signal_time_vect);
        case 'PUPIL_CAMERA'
            [channel_data.pupilFrameOnTimestamps,channel_data.pupilFrameOffTimestamps]=analog_sig_to_timestamps(data(:,params.channels.PUPIL_CAMERA),channel_data.analog_signal_time_vect);
        case 'EEG'
            channel_data.EEG = data(:,params.channels.EEG);
        case 'WATER'
            [channel_data.waterOnTimestamps,channel_data.waterOffTimestamps]=analog_sig_to_timestamps(data(:,params.channels.WATER),channel_data.analog_signal_time_vect);
        case 'LICK'
            [ ~, eventTimes ] = CEDS64ReadEvents( fhand, params.channels.LICK, 100000, 0); % load and convert to our sampling rate
            channel_data.lickTimestamps = double(eventTimes)*dsec;
        otherwise
            warning('Unindetified channel name');
    end
end
end
