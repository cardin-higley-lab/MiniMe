function chain_dual_detrend(id_date_time)
addpath(genpath('/gpfs/ysm/project/cardin/cjb93/MesoProcessing-master'))
addpath('/gpfs/ysm/project/cardin/cjb93/meso-aux-scripts')
BCOSFIREPath = '/gpfs/ysm/project/cardin/cjb93/B-COSFIRE';
addpath(genpath(BCOSFIREPath));
tic

run('../defineIODirs.m'); % input and output directories
run('defineMesoParams.m'); % parameters


currDataDir = fullfile(fixedInputDir,id_date_time);
outputDir = fullfile(fixedOutputDir,id_date_time);
load(fullfile(currDataDir, 'tform_blue.mat'),'tform','R','C');

R=256;
C=256;

disp(['Loading parameters: ' num2str(toc)])



%% load tiffs
tic 
disp(['load tiffs from ' fullfile(currDataDir,'meso') ]);
mov = smart_read_tiffs(fullfile(currDataDir,'meso'),true,params.loading.resizeRatio,-1);

disp(['Loading tiffs: ' num2str(toc)])

%% demixing channels
tic

disp('Demixing channels...')
sigsMov = demix_n_channels(mov, 3);

clear mov;

disp(['demixing channels: ' num2str(toc)])


%% register cameras 
[tform_register,test] = register_image_tform(mean(sigsMov.ch1(:,:,1:100),3),mean(sigsMov.ch3(:,:,1:100),3),BCOSFIREPath);

tform_register_file = fullfile(outputDir,'tform_register.mat');
save(tform_register_file,'tform_register');

% register green frames to blue
sigsMov.ch3 = transform_frames(sigsMov.ch3, tform_register, R, C);

%% align to atlas 
tic

disp('Aligning to template...');
sigsMov.ch1 = transform_frames(sigsMov.ch1, tform, R, C);
sigsMov.ch2 = transform_frames(sigsMov.ch2, tform, R, C);
sigsMov.ch3 = transform_frames(sigsMov.ch3, tform, R, C);



disp(['Aligning to atlas: ' num2str(toc)])

%% detrending
% disp('Detrending...');
tic

% assuming chanel 1 is Neuromodulator signal with potential slow components:
sigsMov.ch1 = reshape(sigsMov.ch1,R*C,[]);
[sigsMov.ch1, ~] = detrend_pixels(sigsMov.ch1, params.detrend);
%sigsMov.ch1 = bsxfun(@minus,bsxfun(@rdivide,sigsMov.ch1,nanmean(sigsMov.ch1, 2)),1);

sigsMov.ch2 = reshape(sigsMov.ch2,R*C,[]);
[sigsMov.ch2, ~] = detrend_pixels(sigsMov.ch2, params.detrend);

sigsMov.ch3 = reshape(sigsMov.ch3,R*C,[]);
[sigsMov.ch3, ~] = detrend_pixels(sigsMov.ch3, params.detrend);



disp(['Detrending: ' num2str(toc)])

%% saving 
tic


sigsMov.ch1 = reshape(sigsMov.ch1,R,C,[]);
sigsMov.ch2 = reshape(sigsMov.ch2,R,C,[]);
sigsMov.ch3 = reshape(sigsMov.ch3,R,C,[]);


mkdir(fullfile(outputDir,'detrended'))
for i = 1:C
    sdet_uv = squeeze(sigsMov.ch1(:,i,:));
    sdet_blue = squeeze(sigsMov.ch2(:,i,:));
    sdet_green = squeeze(sigsMov.ch3(:,i,:));
    save(fullfile(outputDir,'detrended',['detrendedDataUVCol' num2str(i) '.mat']),'sdet_uv');
    save(fullfile(outputDir,'detrended',['detrendedDataBlueCol' num2str(i) '.mat']),'sdet_blue');
    save(fullfile(outputDir,'detrended',['detrendedDataGreenCol' num2str(i) '.mat']),'sdet_green');
end

disp(['Saving: ' num2str(toc)])
