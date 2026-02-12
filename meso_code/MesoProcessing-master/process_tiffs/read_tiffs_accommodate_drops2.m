function [mov, R, C] = read_tiffs_accommodate_drops2(tiffsPath,fs,resizeRatio,varargin)
%%READ_TIFFS_ACCOMMODATE_DROPS load tiff files accounting for dropped
%%frames
% default values for optional arguments:
%   KeepDual: False
%   NumFrames: number of detected tiff files
% usage: mov =
% read_tiffs_accommodate_drops(tiffsPath,offset,frameDiffs,0.5,'KeepDual',true,'NumFrames',1000)
% accommodate wanting to keep both sides of dual channel
% image

    
    
keepDual = false;
if ~isempty(find(strcmp(varargin,'KeepDual'),1))
    keepDual = varargin{find(strcmp(varargin,'KeepDual'),1)+1};
end

if ~isempty(find(strcmp(varargin,'Order'),1))
    Order = varargin{find(strcmp(varargin,'Order'),1)+1};
else
    Order = [1];
end

% get tiffs to load and sort them 
tiffs=dir(fullfile(tiffsPath, '*.tif'));
names = {tiffs.name};
[~,ndx] = natsort(names);
tiffs=tiffs(ndx);

% see what final size of image will be by loading a sample (sniffing)
im = imread(fullfile(tiffs(1).folder,tiffs(1).name));
im = imresize(im, resizeRatio, 'nearest');
[R,C] = size(im);
dual = false;
if R ~= C && ~keepDual
    R=min(size(im));
    C=min(size(im));
    dual = true;
end

n = length(tiffs);
% accommodate loading less frames than total detected
if ~isempty(find(strcmp(varargin,'NumFrames'),1))
    n = varargin{find(strcmp(varargin,'NumFrames'),1)+1};
    if n > length(tiffs)
        % asked for more tiffs than there are. TODO add warning here
        n = length(tiffs);
    end
end

% get timestamp from tiff file
text = fileread(fullfile(tiffs(1).folder,tiffs(1).name));
k = strfind(text,'Time_From_Start = ');
timestamp = text(k+18:k+18+12);
prev_timestamp = str2double(timestamp(1:2))*60*60+str2double(timestamp(4:5))*60+str2double(timestamp(7:8))+str2double(timestamp(9:13));

offset = round(prev_timestamp/(1/fs));




% preallocate movie variable space
% num intended frames = sum(frameDiffs)+offset+1
mov = nan(R,C,sum(frameDiffs)+offset+1,'single');
size(mov)
wb = waitbar(0,'Loading tiffs...');





it = 1+offset;


if dual && ~keepDual
    % iterate through each file loading and performing relevant alterations
    % if dual, we take the side of the image that is the brightest (Current
    % way we do dual imaging)
    for i =1:n
        text = fileread(fullfile(tiffs(i).folder,tiffs(i).name));
        k = strfind(text,'Time_From_Start = ');
        timestamp = text(k+18:k+18+12);
        timestamp = str2double(timestamp(1:2))*60*60+str2double(timestamp(4:5))*60+str2double(timestamp(7:8))+str2double(timestamp(9:13));
        im = imread(fullfile(tiffs(i).folder,tiffs(i).name));
        im = imresize(im, resizeRatio, 'nearest');
%         luminance = zeros(2);
%         luminance(1) = sum(im(:,1:C),'all');
%         luminance(2) = sum(im(:,C+1:end),'all');
%         [~,maxLumInd] = max(luminance);
        maxLumInd = Order(mod(it-1,length(Order))+1);
        im = im(:,((maxLumInd-1)*C+1):(maxLumInd*C));        
        mov(:,:,it) = im;
        if i <= length(frameDiffs)
            it = it + round((timestamp-prev_timestamp)/(1/fs));
            prev_timestamp = timestamp;
        end
        if mod(i,1000) == 0
            disp(['Read ' num2str(i) ' frames']);
        end
        waitbar(i/n,wb,'Loading tiffs...');
    end
else
    % iterate through each file loading and performing relevant alterations
    for i =1:n
        text = fileread(fullfile(tiffs(i).folder,tiffs(i).name));
        k = strfind(text,'Time_From_Start = ');
        timestamp = text(k+18:k+18+12);
        timestamp = str2double(timestamp(1:2))*60*60+str2double(timestamp(4:5))*60+str2double(timestamp(7:8))+str2double(timestamp(9:13));
        im = imread(fullfile(tiffs(i).folder,tiffs(i).name));
        im = imresize(im, resizeRatio);
        mov(:,:,it) = im;
        if i <= length(frameDiffs)
            it = it + round((timestamp-prev_timestamp)/(1/fs));
            prev_timestamp = timestamp;
        end
        waitbar(i/n,wb,'Loading tiffs...');
    end
end