% signal extraction
params.signalsExtraction.firstCaFrame = 1; % is the first frame of tiffs blue
params.signalsExtraction.blueuvRatio = 1; % the ratio between blue and uv imaging frames

% detrending params
params.detrend.filtLen = 1000;
params.detrend.filtcutoff = 0.001;
params.detrend.samplingFreq = 10;
params.detrend.method = 'Butter';
params.detrend.r = 5;
% loading
params.loading.batchSize = 500000000000;
params.loading.batches2load=-1; % load all
params.loading.resizeRatio = 0.5;
R = 256;
C = 256;

% hemo correction
params.hemo.patchSize = 5;
params.hemo.batchSize = 1024;
params.hemo.ch1Blur = 5;
params.hemo.ch2Blur = 5;
params.hemo.ch3Blur = 5;
