function data = reconstructCorrectedData(hemoCorrectedFolder,signalType)
    
    % get all files to read from data dir
    files = dir(fullfile(hemoCorrectedFolder,'*.mat'));
    % find indices that dont correspond to channel we are after
    removeInds = [];
    for i = 1:length(files)
        if ~contains(files(i).name,signalType,'IgnoreCase',true)
            removeInds(end+1) = i;
        end
    end
    
    % remove files not from channel of interest
    files(removeInds) = [];
    
    % peek at data to see how many samples were taken
    peekData = load(fullfile(files(1).folder,files(1).name));
    f = fields(peekData);
    data = nan(256,256,size(peekData.(f{1}),2),'single');
    
    % load data and determine column numbers from file name
    for i = 1:length(files)
        if contains(files(i).name,signalType,'IgnoreCase',true)
            peekData = load(fullfile(files(i).folder,files(i).name));
            fileName = files(i).name;
            startInd = strfind(files(1).name,'Col');
            col = str2double(regexp(fileName(startInd:end),'\d*','Match'));
            data(:,col,:) = single(peekData.(f{1}));
        end
    end

end