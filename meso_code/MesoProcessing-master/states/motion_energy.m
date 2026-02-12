function [motion_energy] = motion_energy(vidDir)

%vidFile = dir(fullfile('K:\rachel\crispr_grin2b\rawData\ro288-01\01_01312024\face\*.avi'));
v = VideoReader(fullfile(vidDir.folder, vidDir.name));

%% how to select roi
% exampleframe = read(v, 1);
% figure();imagesc(exampleframe)
% roi = drawrectangle();
% 
% x = ceil(roi.Position(1));
% y = ceil(roi.Position(2));
% width = ceil(roi.Position(3));
% height = ceil(roi.Position(4));
% testframe = exampleframe(y:y+height,x:x+width,:);
% figure();imagesc(testframe)

%% assess motion energy
tic;
framechange_mean = nan(v.NumFrames-1, 1);
%whiskerchange_mean = nan(v.NumFrames-1, 1);

for i = 1:ceil(v.NumFrames/1000)
    disp([num2str(i) '/' num2str(ceil(v.NumFrames/1000))])

    firstFrame = (i-1)*1000+1;
    lastFrame = i*1000;
    if i == ceil(v.NumFrames/1000)
        lastFrame = v.NumFrames-1;
    end
    
    frames = read(v, [firstFrame lastFrame+1]);
    %whiskers = frames(y:y+height,x:x+width, :, :);

    frame_changes = frames(:,:,:,2:end)-frames(:,:,:,1:end-1);
    framechange_mean(firstFrame:lastFrame) = squeeze(mean(frame_changes, 1:3));
    
    %whisker_changes = whiskers(:,:,:,2:end)-whiskers(:,:,:,1:end-1);
    %whiskerchange_mean(firstFrame:lastFrame) = squeeze(mean(whisker_changes, 1:3));

end
toc;

framechange_mean(isnan(framechange_mean)) = [];
motion_energy = zscore(framechange_mean);

end