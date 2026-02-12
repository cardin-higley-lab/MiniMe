function pre_hemo(id_date_time)

%id_date_time = 'temp\360_12162022';

fprintf([id_date_time ' started processing at %s\n'], datestr(now,'HH:MM:SS.FFF'))

addpath(genpath('../MesoProcessing-master'))

run('../defineIODirs.m'); % input and output directories

outputDir = fullfile(fixedOutputDir,id_date_time);
findslash = strfind(id_date_time,'/');
project = id_date_time(1:findslash(1)-1);
animal = id_date_time(findslash(1)+1:findslash(2)-1);
session = id_date_time(findslash(2)+1:end);


 %% process states if havent..

        if exist(fullfile(outputDir, 'facecam_motion_energy.mat'), 'file')==2
            load(fullfile(outputDir, 'facecam_motion_energy.mat'), 'vid_energy')
        else 
             facemapfile = dir(fullfile(outputDir, '*proc.mat'));
             load(fullfile(outputDir, facemapfile.name), 'proc')
             vid_energy = proc.motSVD{1, 2}(:,1);
        end

        load(fullfile(outputDir, 'final_timestamps.mat'),'spike2_data');
        states = newStates(spike2_data, vid_energy, 2, 3, 10);
        save(fullfile(outputDir, 'state_timestamps.mat'), 'states', '-v7.3');

        statesplot = states_fig(states, spike2_data, vid_energy, animal, session);
        saveas(statesplot, fullfile(outputDir, 'statesplot'));

        behavior = newStates(spike2_data, vid_energy, 2, 0, 0);
        save(fullfile(outputDir, 'behavior_timestamps.mat'), 'behavior', '-v7.3');


%% load raw mats 

disp(['project is: ' project])

if strcmp(project, 'Courtney')
    load(fullfile(outputDir, 'raw_mats.mat'),'gcamp_raw', 'uv_raw');
    blue_raw = gcamp_raw;
elseif strcmp(project, 'CDKL5_dual_spont')
    load(fullfile(outputDir, 'raw_mats.mat'),'rcamp_raw', 'uv_raw', 'blue_raw');
    if exist('rcamp_raw', 'var')
        blue_raw = rcamp_raw;
    end
    disp('success')
else
    load(fullfile(outputDir, 'raw_mats.mat'),'blue_raw', 'uv_raw', 'rcamp_raw');
    if exist('rcamp_raw', 'var')
        blue_raw = rcamp_raw;
    end
 
end

C = 256;


blue_raw = reshape(blue_raw, 256, 256, []);
uv_raw = reshape(uv_raw, 256, 256, []);

% look to see if expression looks ok
temp = blue_raw(:, 14000);
value = prctile(temp(:), 90);
disp(['ninety percent value of raw pixels: ' num2str(value)])
writematrix(value, fullfile(outputDir, 'ninety_percentile_value.txt') )


%% what does the hemo signal look like?

clip = uv_raw(:, :, 1400:1500);

filename = fullfile(outputDir, 'hemo_sig.gif');

for i =1:size(clip,3)
    h=figure('Position', [100 100 900 700]);
    toShow = clip(:,:,i);
    h2 = imagesc(toShow);
    set(h2,'alphadata',~isnan(toShow));
    axis off
    low_lim = quantile(clip(:), 0.1);
    high_lim = quantile(clip(:), 0.99);
    clim([low_lim high_lim]);
    drawnow;
    frame = getframe(h);
    im = frame2im(frame);
    [imind,cm] = rgb2ind(im,256);
    if i == 1
          imwrite(imind,cm,filename, 'gif', 'Loopcount',inf,'DelayTime',.1); %
    else
          imwrite(imind,cm,filename, 'gif','WriteMode','append','DelayTime',.1); %
    end
    close(h)
end



%% Save column-wise
tic


mkNewDir(fullfile(outputDir,'aligned')) 
for batchNum = 1:C/4 
    
    cols = ((colPerBatch*batchNum)-(colPerBatch-1)):(batchNum*colPerBatch);
    sdet_blue = blue_raw(:,cols,:);
    sdet_uv = uv_raw(:,cols,:);

    save(fullfile(outputDir,'aligned',['DataBlueCol' num2str(batchNum) '.mat']),'sdet_blue');
    save(fullfile(outputDir,'aligned',['DataUVCol' num2str(batchNum) '.mat']),'sdet_uv');
end

disp('finished saving')
clear blue_aligned uv_aligned

disp(['Saving: ' num2str(toc)])

end
