%% big iteration through data

% adjusted from statePartialCorrelation.m 
% editted to correlate for loco state ONLY, no faceMap data 

clearvars
load('parcel_Obj.mat');

% spike2data in normal location
masterFolder = 'C:\Users\rachel\Desktop\PreprocessedData\Ketamine';
subFolders = dir([masterFolder '\*\*']);

partialCorrData = struct;
for folder = 1:length(subFolders)
    
    % for testing 
    folder='pv1_day1';
    final_timestampsObj = load('C:\Users\rachel\Desktop\meso_rlo\fixedTimestamps\pv1_day1\final_timestamps.mat')
    load('C:\Users\rachel\Desktop\PreprocessedData\Ketamine\PV1\pv1_day1\final_dFoF_parcels.mat')
    
    % ignore non data folders (.,..)
    if length(subFolders(folder).name)<3 
        continue;                                         
    end

    % ignore subfolders that didnt process correctly 
    folders2ignore = ["pv2_day2", "syn1_day1", "pv1_day4"];
    if contains(subFolders(folder).name, folders2ignore) == 1
        continue;
    end
    
    disp(['Processing: ' subFolders(folder).name]);
    
    % get path to load data
    dataFolderPath = [subFolders(folder).folder '\' subFolders(folder).name];
    
    % load parcels dff and timestamp data
    load(fullfile(dataFolderPath,'final_dFoF_parcels.mat'));
    
    % load appropriate timestamp data 
    final_timestampsObj = load(['C:\Users\rachel\Desktop\meso_rlo\fixedTimestamps\' subFolders(folder).name '\final_timestamps']); 
    spike2_data = final_timestampsObj.spike2_data; % extract spike2data 
   
    %get timestamps for states
    states = locoTimestamps(spike2_data);
   
    %get the number of parcels and the number of frames 
    [numParcels,numSamples] = size(dFoF_parcells.blue);
    
    %take midpoint of exposure time as timepoint
    mesoTimestamps = spike2_data.mesoFrameOnTimestamps(250:2:end-250);
    mesoTimestamps = mesoTimestamps(1:length(dFoF_parcells.blue));
    
    %find parcels w NaNs (e.g. colliculi) and remove them
    nanlessParcelInds = ~isnan(dFoF_parcells.blue(:,1));
    nanlessParcels = dFoF_parcells.blue(nanlessParcelInds,:);
   
 
    % weights the parcels according to size??? finding the mean (but not by parcel..?
    temp1 = sum(parcel_obj.indicators(:,:,nanlessParcelInds),[1 2]); % returns 1x1x52 double
    
    
    %   S = SUM(X,VECDIM) operates on the dimensions specified in the vector VECDIM. 
    %   SUM(X,[1 2]) operates on the elements contained in the first and second dimensions of X.
    
    temp2 = sum(parcel_obj.indicators(:,:,nanlessParcelInds),'all'); % = 34253
    
    %   S = SUM(X,'all') sums all elements of X.
    
    
    weightings = squeeze(temp1/temp2);
    meanTrace = sum(nanlessParcels.*weightings,1); % not sure why you would want a single trace..?
    
    %
    stateFields = fields(states);
    numStates = length(stateFields)/2;
    
    
    % for locomotion
    currOnTimestamps = states.locoOn;
    currOffTimestamps = states.locoOff;
    
    %get frames that fall within the bounds of the current state
    %stateInds = zeros(numSamples,1); % this doesnt work.. 
    for ii = 1:length(currOffTimestamps)
        stateInds(mesoTimestamps >= currOnTimestamps(ii) & mesoTimestamps<=currOffTimestamps(ii)) = true;
    end
    
    %calculate partial correlation and save in standard format   
    partialCorrMats = nan(numParcels,numParcels,numStates); %mat to hold our partial correlation data
    
    
    % variables to run partial correlation on
    temp3 = nanlessParcels(:,stateInds)'; %4973x52 double --> for each parcel, what the fluorescent signal was for each locoOn frame 
    temp4 = meanTrace(stateInds)'; % 4973x1 double --> mean fluorescence for each locoOn frame 
    
    % correlations --> why correlating to mean trace ???
    % ok so (se below) need to control for mean trace (weightings) for some
    % reason but not sure why..? 
    partialCorrMats(nanlessParcelInds,nanlessParcelInds,1) = partialcorr(temp3 , temp4);
    
    
    %   RHO = PARTIALCORR(X,Z) returns the sample linear partial correlation
    %   coefficients between pairs of variables in X, controlling for the
    %   variables in Z.  X is an N-by-P matrix, and Z an N-by-Q matrix, with rows
    %   corresponding to observations, and columns corresponding to variables. RHO
    %   is a symmetric P-by-P matrix. 
    
   
    
    session = subFolders(folder).name;
    partialCorrData.(session) = struct;
    partialCorrData.(session).corrs = partialCorrMats;
    partialCorrData.(session).session = session;
end


