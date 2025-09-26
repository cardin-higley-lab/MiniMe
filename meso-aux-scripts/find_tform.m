function post_hemo_regular(mouseID, date)
%% load id_date_time specific info
tic

addpath(genpath('../MesoProcessing-master'))
addpath('../meso-aux-scripts')

run('defineMesoParams.m'); % parameters

main_folder = '/vast/palmer/scratch/higley/hd362/HD_Mouse_Training';

mouseID = 'mouse2';
date = '0211';

mouseID(1) = upper(mouseID(1));

outputFolder = dir(fullfile(main_folder,strcat('*', mouseID, '_', date, '*'), '*', '*Miniscope'));
outputDir = outputFolder.folder;


%load(fullfile(currDataDir, 'tform_blue.mat'),'R','C');
R=256;
C=256;

disp(['Loading parameters: ' num2str(toc)])

%% Load column data into a full data matrix
tic
% define the folder where data is kept
dataDir = outputDir;

% load allen parcels
load('/gpfs/gibbs/project/higley/hd362/MesoProcessing-master/allen_parcel/parcels_updated12522.mat', 'allen_parcels');

% load hemo corrected data
frameRange = [10000, 13000];
data = loadColumnwise(fullfile(dataDir,'hemoCorrectedSig'),'blue', frameRange);

toc

%% Find the tform and save the matrix
movVar = var(data,[],3,'omitnan');
lowPrctile = prctile(movVar(:),5);
highPrctile = prctile(movVar(:),95);
movVar(movVar<lowPrctile) = lowPrctile;
movVar(movVar>highPrctile) = highPrctile;
reactiveAllenAlignment(movVar)

tformFile = fullfile(main_folder, 'process', strcat('tform', mouseID, '.mat'));
save(tformFile, 'tform');

[mouseFolder,~,~] = fileparts(dataDir);
tformFile = fullfile(mouseFolder, strcat('tform', mouseID, '.mat'));
save(tformFile, 'tform');

%% Load tform.mat and transform the data (old version)
tformFile = fullfile(main_folder, 'process', strcat('tform', mouseID, '.mat'));

load(tformFile);

data_tform = transform_frames(data,tform);

%% Load tform.mat and tranform the data (full version)

% tformFolder = dir(fullfile(main_folder, strcat('*', mouseID, '_', date, '*'), 'tform*'));
% tformFile = fullfile(tformFolder(1).folder, tformFolder(1).name);
% 
% load(tformFile);
% 
% tform = simtform2d(tform.Scale, tform.RotationAngle, tform.Translation);
% 
% data_tform = transform_frames(data,tform);

%% Use Mask to block out non brain pixels

mask = allen_parcels.CombinedParcells>0;

% mask points that don't fall on the brain and zscore those that do
for i =1:256
   data_tform(~mask(:,i),i,:) = nan;
   % data_tform(mask(:,i),i,:) = normalize(data_tform(mask(:,i),i,:),3);
end
%% Parcellate the data

parcellated = parcels_by_Allen_atlas(reshape(data_tform,256^2,[]),allen_parcels);

parcelFile = strcat('parcellated_', mouseID, '_', date, '.mat');
% save parcellated data
save(fullfile(dataDir, parcelFile),'parcellated')

%% Save the data_tform into mat

data_tformFile = strcat('Neural_', mouseID, '_', date, '.mat');
data_tformfilePath = fullfile(dataDir, data_tformFile);

% Check if the file already exists
if ~exist(data_tformfilePath, 'file')
    % If it does not exist, save the data
    save(data_tformfilePath, 'data_tform', '-v7.3');
    disp(['File saved: ', data_tformfilePath]);
else
    disp(['File already exists: ', data_tformfilePath]);
end

%% write data into avi video
% Sample data (replace with your actual data)
data = data_tform;
% data = clips_mean;
% Define the range of pixel values you want to display
min_val = -3;  % Minimum pixel value to display
max_val = 3;  % Maximum pixel value to display

video_filename = fullfile(dataDir, strcat('Neural_HC_R', num2str(max_val), '_', mouseID, '_', date));

%video_filename = 'W:\Hao\Globus\2024_11_04_15_29\clip_ave_v2'
% Create a VideoWriter object, specify the filename and frame rate
v = VideoWriter(video_filename);  % Specify the file name and format
v.FrameRate = 30;  % Set the frame rate to 30 Hz

% Open the video file for writing
open(v);

% Loop through each frame, adjust pixel values, and write to video
for i = 1:size(data, 3)
    frame = data(:, :, i);  % Extract the ith frame

    % Normalize the frame to the range [0, 1] based on the specified range
    adjusted_frame = (frame - min_val) / (max_val - min_val);

    % Clip values outside the range [0, 1]
    adjusted_frame(adjusted_frame < 0) = 0;
    adjusted_frame(adjusted_frame > 1) = 1;

    % Convert the adjusted frame to RGB using a colormap
    rgb_frame = ind2rgb(uint8(adjusted_frame * 255), colormap('parula'));

    % Write the RGB frame to the video
    writeVideo(v, rgb_frame);

end

% Close the video file
close(v);

disp('Video saved successfully.');

end
