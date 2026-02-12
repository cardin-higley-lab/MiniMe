%% user-selected input parameters 
params.svdFact = 14;
params.TmpFile = 'tmp.tif';
params.cedpath = 'D:\meso_demo\MesoProcessing-master\process_spike2\CEDS64ML';
params.Loc=location; %recording location, 'B449' or 'F238'
params.fsspike2=5000;
params.fsimaging=10; 
params.pupilSR=10; 
params.batchSize = 3000;
params.batches2load=10000;
params.regressUV=1; %option of regressing every pixel by UV 
params.moveLocal=0;%if you want to move the tiff file to a local drive before reading 
params.drawPlots=0;%plot all event channels and imaging channel, for test purposes
params.makeMovies=0;%make tif df/f movies for a segement of the session 
params.imageFormat='cxd';%choose between cxd or tif 
% visusal stim 
params.visStimAn=1; % 1 if vis stim are presented, 0 if not 
if strcmp(sessid,'Vis'), params.visStimAn=1;end 
params.airpuffAn=0; % if airpuffs or electrical stim are presented, indicate 1; 
if strcmp(sessid,'Spont'), params.airpuffAn=1;end 
params.visDur=2;
params.visITI=5; 
%grabs/rcamp vs grabs only
params.signalsExtraction.firstCaFrame = 1;
params.signalsExtraction.sigs = color;% 'blueuv' or 'RCaMP_AC'
params.signalsExtraction.blueuvRatio = 1;
params.resizeRatio = 0.5;
% detrending and spatial filtering params
params.deterend.filtLen = 100;
params.deterend.filtcutoff = 0.01;
params.deterend.method = 'FIR';
params.spatiotemporalfilter=1; %set it to 1 if you want to do spatial filtering  otherwise 0 

%run params
params.minRunDuration=5; %in seconds
params.minSitDuration=10; %in seconds 
params.ITITime=5; %  dead time after events in seconds  

%analysis window for events (airpuff and visstim)
params.preEventWin=2;
params.postEventWin=5;
params.baselineWin=2; %use this period to calculate baseline during within trial normalization 
%% choose right parameters depending on the recording location 
if strcmp(params.Loc,'B449')
    params.angle2rotate=-180;
    % spike2 channels for bcmm 
    channels.BLUE = 1;%blue led 
    channels.UV = 2;%uv led 
    channels.FRAMETICKS = 8;%green/uv mesocam triggers
    channels.PHOTO_DIODE = 4;%visual stim trigger
    channels.WHEEL = 5;
    channels.AIR_PUFF = 6; %this is either airpuff channel or electrical stim channel 
    channels.PUPIL_CAMERA = 7;
    channels.RED_MESO = 3;%red mesocam trigger
    channels.GREEN=9;%green LED
    channels.EEG = 10;%eeg continous signal 
elseif strcmp(params.Loc,'F238')
    params.angle2rotate=90;
    % spike2 channels for cardin 238
    channels.BLUE = 1;%blue led 
    channels.UV = 2;%uv led 
    channels.FRAMETICKS = 3;%green/uv mesocam triggers
    channels.PHOTO_DIODE = 4;%visual stim trigger
    channels.WHEEL = 5;
    channels.AIR_PUFF = 6; %this is either airpuff channel or electrical stim channel 
    channels.PUPIL_CAMERA = 7;
    channels.RED_MESO = 3;%red mesocam trigger
    channels.WATER = 8;% water delivery, analog
    channels.LICK = 9; % lick events
end 
params.channels=channels;
