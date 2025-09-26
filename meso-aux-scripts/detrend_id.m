function detrend_id(id_date_time, batchNum)
addpath(genpath('../MesoProcessing-master'))
addpath(genpath('../meso-aux-scripts'))
batchNum = str2double(batchNum);

%disp(['Argument 1: ', mouseID])
%disp(['Argument 2: ', date])

tic
run('defineMesoParams.m'); % parameters

main_folder = '/vast/palmer/scratch/higley/hd362/HD_Mouse_Training';

% id_date_time = 'SP_MouseF_0820_pre';

outputFolder = dir(fullfile(main_folder,id_date_time, '*', 'My_V4_Miniscope'));
outputFolder = outputFolder.folder;
outputDir= fileparts(outputFolder);

%load(fullfile(currDataDir, 'tform_blue.mat'),'tform','R','C');

R = 256;
C = 256;
n=2;

disp(['Loading parameters: ' num2str(toc)])

tic
for ch = 1:n
    load(fullfile(outputDir, 'RawDemixed',['RawDemixedCH' num2str(ch) 'Col' num2str(batchNum) '.mat']),'column');
    sigsMov.(['ch' num2str(ch)]) = column;
end

disp(['loading channels: ' num2str(toc)])

%% detrending channels
tic 
for ch = 1:n
    sigsMov.(['ch' num2str(ch)]) = detrend_pixels(sigsMov.(['ch' num2str(ch)]), params.detrend);
end

disp(['Detrending: ' num2str(toc)])

%% save
tic

% Check if 'detrended' folder exists, if not, create it
if ~isfolder(fullfile(outputDir, 'detrended'))
    mkdir(fullfile(outputDir, 'detrended'));
end

% Loop through each channel and save the data
for ch = 1:n
    column = sigsMov.(['ch' num2str(ch)]);
    
    % Save each column of data in -v7.3 format
    save(fullfile(outputDir, 'detrended', ['detrendedDataCH' num2str(ch) 'Col' num2str(batchNum) '.mat']), 'column', '-v7.3');
end

disp(['Saving: ' num2str(toc)])
%saveColumnwise(sigsMov.ch1,fullfile(outputDir,'detrended'),'detrendedDataCH1')
%saveColumnwise(sigsMov.ch2,fullfile(outputDir,'detrended'),'detrendedDataCH2')


