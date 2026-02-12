
load('D:\meso_demo\CDKL5\ro076-123\123_04292022\123_04292022-04292022160749_proc.mat')
load('D:\meso_demo\CDKL5\ro076-123\123_04292022\final_timestamps.mat')


facemapLength = length(proc.motSVD{1, 2}(:,1));
pupilLength = length(spike2_data.pupilFrameOffTimestamps);

if  pupilLength - facemapLength > 28 % 28 bc always drops last 27 and 1 extra for slack
    error('significant frame drops exist in facemap data')
end

% check for low face bouts; error if none and normalizing by low face
pc = proc.motSVD{1, 2}(:,1);
if skewness(pc) < 0.5
    pc = -pc;
end
face_data = pc;

hands = proc.motSVD{1, 3}(:,1);
if skewness(hands) < 0.5
    hands = -hands;
end



states = stateTimestamps_grooming(spike2_data, pc, hands);

statesplot = states_fig(states, spike2_data, pc, hands);

sweyta = sweyta_states(spike2_data, pc, hands);

sweyta_statesplot = states_fig(sweyta, spike2_data, pc, hands);



%% removed from stateTimestamps
% boolResult = timestampsIntersect(spike2_data.pupilFrameOnTimestamps,spike2_data.pupilFrameOffTimestamps, sitOn_final, sitOff_final);
% boolResult = boolResult(1:length(pc));
% pc_bool = normalize(pc);
% pc_bool(boolResult == 0) = NaN;
% nonNAs=(pc_bool(~isnan(pc_bool)));