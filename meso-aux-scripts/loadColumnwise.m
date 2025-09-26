function data = loadColumnwise(hemoCorrectedFolder, signalType, frameRange)
    % Load 2D data from .mat files, with row as pixel number and column as frames.
    % Only loads frames within the specified frame range using matfile for efficiency.
    %
    % Inputs:
    %   hemoCorrectedFolder - Path to the folder containing the .mat files.
    %   signalType - The type of signal to filter the files by (e.g., 'EEG', 'EOG').
    %   frameRange - A 2-element vector specifying the start and end frame indices (e.g., [start_frame, end_frame]).
    %
    % Output:
    %   data - A 3D matrix of type 'single' (256x256x(num_frames_in_range)) containing the selected data.

    % Get all files in the directory
    files = dir(fullfile(hemoCorrectedFolder, '*.mat'));
    
    % Find files that do not match the signal type in their name
    removeInds = [];
    for i = 1:length(files)
        if ~contains(files(i).name, signalType, 'IgnoreCase', true)
            removeInds(end+1) = i;
        end
    end
    % Remove files that do not match the signal type
    files(removeInds) = [];
    
    % Peek at data from the first file to determine the number of frames
    peekData = load(fullfile(files(1).folder, files(1).name));
    f = fields(peekData);
    
    % Assuming data is in the first field, and it has the form [256, 256, num_frames]
    % We'll extract the number of frames from the size of the data in the second dimension
    numFrames = size(peekData.(f{1}), 2);
    
    % If frameRange is provided, we check that it's valid
    if nargin < 3 || isempty(frameRange)
        % If no frame range is specified, use the entire set of frames
        frameRange = [1, numFrames];
    else
        % Ensure frameRange is within the bounds
        if frameRange(1) < 1 || frameRange(2) > numFrames || frameRange(1) > frameRange(2)
            error('Invalid frame range specified.');
        end
    end
    
    % Initialize the data matrix: 256x256 for the image dimensions, and frames in the range specified
    data = nan(256, 256, diff(frameRange) + 1, 'single');
    
    % Load the data and fill the matrix using matfile for partial loading
    for i = 1:length(files)
        if contains(files(i).name, signalType, 'IgnoreCase', true)
            % Create a matfile object to load only the required frames
            matObj = matfile(fullfile(files(i).folder, files(i).name));
            
            % Find the column number from the file name (assuming 'ColX' in the name)
            startInd = strfind(files(i).name, 'Col');
            col = str2double(regexp(files(i).name(startInd:end), '\d*', 'Match'));
            
            % Load only the frames in the specified range (rows x selected frames)
            signalData = matObj.(f{1})(:, frameRange(1):frameRange(2));
            
            % Assign the data to the appropriate column in the 3D matrix
            data(:, col, :) = single(signalData);
        end
    end
end
