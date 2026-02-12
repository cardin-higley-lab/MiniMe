load('D:\meso_demo\rawData\ro080-22730\22730_03192022\face\22730_face-03192022154547_proc.mat')
load('D:\meso_demo\CHD8\ro080-22730\22730_03192022\final_timestamps.mat')

pc = proc.motSVD{1, 2}(:,1);
if skewness(pc) < 0.5
    pc = -pc;
end

hands = proc.motSVD{1, 3}(:,1);
if skewness(hands) < 0.5
    hands = -hands;
end

handmask1 = proc.uMotMask{1,3}(:,:,1);
handmask1_image = figure(); imagesc(handmask1);

% sampling_freq = round(1/median(diff(spike2_data.pupilFrameOnTimestamps)));
% [b,a] = butter(2, 0.1/(sampling_freq/2)); % 0.1 = cutoff frequency
% pc = filtfilt(b,a,double(pc));

states = stateTimestamps_grooming(spike2_data, pc, hands);


