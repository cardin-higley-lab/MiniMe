function hemo_correct_id(id_date_time, batchNum)

batchNum = str2double(batchNum);

%% load id_date_time specific info
tic

addpath(genpath('../MesoProcessing-master'))
addpath('../meso-aux-scripts')

run('defineMesoParams.m'); % parameters

main_folder = '/vast/palmer/scratch/higley/hd362/HD_Mouse_Training';

% id_date_time = 'SP_MouseF_0820_pre';

outputFolder = dir(fullfile(main_folder,id_date_time, '*', 'My_V4_Miniscope'));
outputFolder = outputFolder.folder;
outputDir= fileparts(outputFolder);

%load(fullfile(currDataDir, 'tform_blue.mat'),'R','C');
R=256;
C=256;

disp(['Loading parameters: ' num2str(toc)])


%% load detrended data
% only load data we will use. 
% assumes patchSize variable is really half - 1 actual patch size
tic

% columns per batch based on batch size
colPerBatch = ceil(params.hemo.batchSize/C);

% specific columns to process and load, and then corresponding pixels of
% our reduced dataset (compensates for being on edge
colToProcess = (batchNum-1)*colPerBatch+1:batchNum*colPerBatch;
colToLoad = (batchNum-1)*colPerBatch-params.hemo.patchSize+2:(batchNum)*colPerBatch-1+params.hemo.patchSize;
colToLoad(colToLoad<1) = [];
colToLoad(colToLoad>C) = [];
pixelsToProcess = (find(colToLoad == colToProcess(1))-1)*C+1:(find(colToLoad == colToProcess(end)))*C;



% load one column to check how many timepoints there are
load(fullfile(outputDir,'detrended',['detrendedDataCH1Col' num2str(colToLoad(1)) '.mat']),'column');
ch1 = column;
% preallocate space for loading data
det_ch1 = zeros(R,length(colToLoad),size(ch1,2));
det_ch2 = zeros(R,length(colToLoad),size(ch1,2));
% read uv and blue data from each column and add to data matrix
for i =1:length(colToLoad)
    load(fullfile(outputDir,'detrended',['detrendedDataCH1Col' num2str(colToLoad(i)) '.mat']),'column');
    assert(size(column,2) == size(det_ch1,3), ...
    sprintf('Column %d has inconsistent timepoints: expected %d, got %d', ...
    i, size(det_ch1,3), size(column,2)));

    det_ch1(:,i,1:length(column)) = column;

    load(fullfile(outputDir,'detrended',['detrendedDataCH2Col' num2str(colToLoad(i)) '.mat']),'column');
    assert(size(column,2) == size(det_ch1,3), ...
    sprintf('Column %d has inconsistent timepoints: expected %d, got %d', ...
    i, size(det_ch1,3), size(column,2)));

    det_ch2(:,i,1:length(column)) = column;
end

disp(['Loading data: ' num2str(toc)])


%% Remove Hemodynamics
tic

% reshape data to match what regress_filter wants
det_ch1 = reshape(det_ch1,R*length(colToLoad),[]);
det_ch2 = reshape(det_ch2,R*length(colToLoad),[]);

% fit_brain_mask_file = fullfile(outputDir,'fit_brain_mask.mat');
% load(fit_brain_mask_file,'fit_brain_mask');
% mask middle bar

load('parcells_updated121519.mat','parcells_new');
brain_mask = true(R,C);
%brain_mask(:,256/2-6:256/2+6) = 0;


% only pass part of mask that corresponds to pixels in this batch
brain_mask = brain_mask(:,colToLoad);

flat_mask = reshape(brain_mask,1,[]);

dropped_frame_size = size(det_ch1(flat_mask,:),1);


%% regress first channel
nanlessMask = sum(isnan(det_ch2(flat_mask,:)),1) ~= dropped_frame_size &  sum(isnan(det_ch1(flat_mask,:)),1) ~= dropped_frame_size;
disp(size(det_ch1(flat_mask,nanlessMask)))

disp(size(det_ch2(flat_mask,nanlessMask)))

[std_blue, std_uvblue, blue_reg]  = regress_filter(...
    det_ch1(flat_mask,nanlessMask), ...
    det_ch2(flat_mask,nanlessMask), ...
    brain_mask, params.hemo.patchSize, pixelsToProcess);

nblue_reg = nan(size(blue_reg,1),size(det_ch1,2));
nblue_reg(:,nanlessMask) = blue_reg;
blue_reg = nblue_reg;

disp(['Hemodynamic correction: ' num2str(toc)])

%% Save
tic

blue_reg = reshape(blue_reg, R, colPerBatch, []);
std_blue = reshape(std_blue, R, colPerBatch);
std_uvblue = reshape(std_uvblue, R, colPerBatch);

% Create directories if they don't exist
mkdir(fullfile(outputDir, 'hemoCorrectedSig'))
mkdir(fullfile(outputDir, 'hemoCorrectedNoise'))

% Loop through columns
for i = 1:colPerBatch
    % Extract the relevant slices of the data
    sblue_reg = squeeze(blue_reg(:, i, :));
    sstd_uvblue = squeeze(std_uvblue(:, i));
    sstd_blue = squeeze(std_blue(:, i));
    
    % Save in -v7.3 format for partial loading support
    save(fullfile(outputDir, 'hemoCorrectedSig', ['hemoCorrectedSigBlueCol' num2str(colToProcess(i)) '.mat']), 'sblue_reg', '-v7.3');
    save(fullfile(outputDir, 'hemoCorrectedNoise', ['noiseUVBlueCol' num2str(colToProcess(i)) '.mat']), 'sstd_uvblue', '-v7.3');
    save(fullfile(outputDir, 'hemoCorrectedNoise', ['noiseBlue' num2str(colToProcess(i)) '.mat']), 'sstd_blue', '-v7.3');
end

disp(['Saving: ' num2str(toc)])

