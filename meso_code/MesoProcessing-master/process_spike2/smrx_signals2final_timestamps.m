%need to upload allen parcellation traces, final_timestamps, and
%smrx_signals
% for converting smrx_signals.mat (sweyta/andrew old code) to
% final_timetamps structure

project = 'ACh';
animal = 'ach_AM10';
session = '11232019_grabAM10_spont';

masterDir = 'C:\Users\CardinLab\Desktop\meso_demo\';
saveoutputpath = fullfile(masterDir, project, animal, session);


%% need to convert :
spike2_data.blueOnTimestamps = timing.bluestart;
spike2_data.blueOffTimestamps = timing.blueend;
spike2_data.uvOnTimestamps = timing.uvstart;
spike2_data.uvOffTimestamps = timing.uvend;
spike2_data.mesoFrameOnTimestamps = timing.mesostart;
spike2_data.mesoFrameOffTimestamps = timing.mesoend;
spike2_data.wheelOn = timing.wheelOn;
spike2_data.wheelOff = timing.wheelOff;
spike2_data.pupilFrameOnTimestamps = timing.pupilcamstart;
spike2_data.pupilFrameOffTimestamps = timing.pupilcamend;


%% sanity check 

figure();
plot(spike2_data.analog_signal_time_vect,spike2_data.wheelSpeed); 
ylabel('wheel speed');
xline(spike2_data.wheelOn, 'g');
xline(spike2_data.wheelOff, 'r');

%% dff @ locoOnset 
load('parcells_updated121519.mat','parcells_new');

fs = 10; % meso sampling rate (fs = frames/second)
preSeconds = 2; % number of seconds of data we want before the event
postSeconds = 2; % number of seconds of data we want after the event

locoOnset = figure; hold on;
for i = 1:length(1:2:56)
    parcel_names = parcells_new.description(1:2:56);
    colors = distinguishable_colors(length(parcel_names));
    t = -preSeconds:(preSeconds+postSeconds)/(preSeconds*fs+postSeconds*fs):postSeconds; % get time vect for plotting
    %spk2 = spike2_data;
    parcelMean = parcels4locoOnset(parcels_time_trace, spike2_data.blueOnTimestamps, spike2_data.blueOffTimestamps, spike2_data.wheelOn, fs, preSeconds, postSeconds);
    plot(t,parcelMean(i,:), 'Color' , colors(i,:));
end
ylabel('dFoF (%)');
xlabel('Time (s)');
xline(0);
legend(parcel_names, 'Location', 'northwest');

saveas(locoOnset, fullfile(saveoutputpath, 'dffLocoOnset'));

%%
save(fullfile(saveoutputpath, 'final_timestamps.mat'), 'spike2_data', '-v7.3');

