function hemo_correct(id_date_time, batchNum)

tic

% id_date_time = 'WH_MouseK_0703_WH01_3';
% batchNum = '32';

addpath(genpath('../MesoProcessing-master'))
run('../defineIODirs.m'); % input and output directories
%patchSize = 9;colPerBatch = 4;delay_length = 500;
outputFolder = dir(fullfile(fixedOutputDir,id_date_time, '**', 'aligned'));
outputFolder = outputFolder.folder;
outputDir= fileparts(outputFolder);


R=256;C=256;
batchNum = str2double(batchNum);

disp(['Loading parameters: ' num2str(toc)])

%% load detrended data
tic

% specific cols to process... corresponding pixels of reduced dataset (compensates for being on edge)
colToProcess = (batchNum-1)*colPerBatch+1:batchNum*colPerBatch;
colToLoad = (batchNum*colPerBatch-(colPerBatch-1))-(patchSize-1)/2:(batchNum)*colPerBatch+(patchSize-1)/2;
colToLoad(colToLoad<1) = [];
colToLoad(colToLoad>C) = [];

halfPatch = (patchSize-1)/2;
numColsNeeded = ceil(halfPatch/colPerBatch); 

batchesToLoad = (batchNum-numColsNeeded):(batchNum+numColsNeeded);
batchesToLoad(batchesToLoad<1) = [];
batchesToLoad(batchesToLoad>(C/colPerBatch)) = [];

% load one to get size and preallocate
uvdata = load(fullfile(outputDir,'aligned',['DataUVCol' num2str(batchesToLoad(1)) '.mat']),'column');
sdet_uv = permute(uvdata.column, [2, 1, 3]);
alignedUVColData = nan(size(sdet_uv, 1), size(sdet_uv, 2)*length(batchesToLoad), size(sdet_uv,3));
bldata = load(fullfile(outputDir,'aligned',['DataBlueCol' num2str(batchesToLoad(1)) '.mat']),'column');
sdet_blue = permute(bldata.column, [2, 1, 3]);
alignedBLColData = nan(size(sdet_blue, 1), size(sdet_blue, 2)*length(batchesToLoad), size(sdet_blue,3));

for i = 1:length(batchesToLoad)
    cols = ((colPerBatch*i)-(colPerBatch-1)):(i*colPerBatch);
    uvdata = load(fullfile(outputDir,'aligned',['DataUVCol' num2str(batchesToLoad(i)) '.mat']),'column');
    sdet_uv = permute(uvdata.column, [2, 1, 3]);
    alignedUVColData(:,cols,:) = sdet_uv;
    bldata = load(fullfile(outputDir,'aligned',['DataBlueCol' num2str(batchesToLoad(1)) '.mat']),'column');
    sdet_blue = permute(bldata.column, [2, 1, 3]);
    alignedBLColData(:,cols,:) = sdet_blue;
end

pixelsToProcess = (find(colToLoad == colToProcess(1))-1)*C+1:(find(colToLoad == colToProcess(end)))*C;

minLen = min(size(alignedUVColData,3),size(alignedBLColData,3));
det_blue = alignedBLColData(:,:,1:minLen);
det_uv = alignedUVColData(:,:,1:minLen);


disp(['Loading data: ' num2str(toc)])

%% threshold for bad pixels, medfilt hemo signal, and mean subtract

% mean subtract signal
signal_raw = cat(3, nan(size(det_blue(:, :, 1:delay_length))), det_blue(:, :, delay_length+1:end)); % converts beginning of session to nans
signal_mean = mean(signal_raw, 3, 'omitnan'); % finds the pixelwise-mean of the matrix
signal_sub = (bsxfun(@minus, signal_raw, signal_mean)); % subtracts the mean from detrended mat
signal_z = reshape(signal_sub, size(signal_sub,1)*size(signal_sub,2), []);


% process hemo signal and then mean subtract
hemo_raw = cat(3, nan(size(det_uv(:, :, 1:delay_length))), det_uv(:, :, delay_length+1:end)); % converts beginning of session to nans
hemo_mean = mean(hemo_raw, 3, 'omitnan'); % finds the mean of the detrended matrix

% flag over or undersaturated pixels and set as nan
hemo_raw_flat = reshape(hemo_raw, size(hemo_raw, 1) * size(hemo_raw, 2), []);
too_large = quantile(hemo_mean(:), 0.95) < hemo_mean;
too_small = quantile(hemo_mean(:), 0.05) > hemo_mean;
hemo_raw_flat(too_small, :) = NaN;
hemo_raw_flat(too_large, :) = NaN;

% use linear interpolation to fill in nan values
temp3 = fillmissing(hemo_raw_flat, 'linear');
temp3 = reshape(temp3, size(hemo_raw, 1), size(hemo_raw, 2), []);

% apply small median filter 
temp4 = nan(size(temp3));
for i = 1:size(temp4, 3)
    temp4(:,:,i) = medfilt2(temp3(:,:,i), [3 3]);
end

hemo_mean = mean(temp4, 3, 'omitnan'); % finds the mean of the detrended matrix
hemo_sub = (bsxfun(@minus, temp4, hemo_mean)); % subtracts the mean from detrended mat
hemo_z = reshape(hemo_sub, size(hemo_sub,1)*size(hemo_sub,2), []);


%% Remove Hemodynamics
tic

%sub_brain_mask = brain_mask(:,colToLoad);
%flat_mask = reshape(sub_brain_mask,1,[]);

% run hemo
tic
[~, ~, blue_reg]  = regress_filter(...
    signal_z(:, delay_length+1:end), ... 
    hemo_z(:, delay_length+1:end), ... 
    ones(256, 12), patchSize, pixelsToProcess);
toc

output_sig = cat(2, nan(size(blue_reg, 1), delay_length), blue_reg); 


disp(['Hemodynamic correction: ' num2str(toc)])

%% Undo Z-scoring

output = reshape(output_sig,R,colPerBatch,[]); %implay(output(:,:,500:end))

cols_wanted = nan(1, length(colToProcess));
for col = 1:length(colToProcess)
    cols_wanted(col) = find(colToProcess(col)==colToLoad);
end

blue_mean_crop = signal_mean(:, cols_wanted);
new_output = (bsxfun(@plus, output, blue_mean_crop)); 


%% Save
tic

sblue_reg = new_output;
mkNewDir(fullfile(outputDir,'hemoCorrectedSigDff'))
save(fullfile(outputDir,'hemoCorrectedSigDff',['hemoCorrectedSigCol' num2str(batchNum) '.mat']),'sblue_reg');

% sstd_blue = reshape(std_blue,R,colPerBatch);
% sstd_uv = reshape(std_uv,R,colPerBatch);
% mkNewDir(fullfile(outputDir,'hemoCorrectedNoise'))
% save(fullfile(outputDir,'hemoCorrectedNoise',['noise' num2str(batchNum) '.mat']),'sstd_uv','sstd_blue');


disp(['Saving: ' num2str(toc)])

