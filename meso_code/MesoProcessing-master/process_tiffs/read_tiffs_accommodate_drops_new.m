function [mov, R, C, dual] = read_tiffs_accommodate_drops_new(tiffsPath,fs,ResizeRatio, Order, numFrames)
%%READ_TIFFS_ACCOMMODATE_DROPS load tiff files froma a directory accounting for dropped frames
%
% Required input arguments:
%    tiffsPath: The directory where the tiff files can be found. 
%    fs: The sampling frequency of imaging (overall fs, not per channel)
%    ResizeRatio: 0.5 is standard (512 to 256)
%    Order: [1] 
%    numFrames: number of frames to read in, -1 if all 
%%

%tiffsPath = 'C:\Users\CardinLab\Desktop\rawData\ro168-474\test\meso';
%fs = 30;
%ResizeRatio=0.5;
%Order=[1 1 2];
%numFrames=-1;

%%
% quick check to make sure that the order and fs makes sense 
if length(Order)*10 ~= fs
    error('Order length and frame rate disagree')
end

% get tiffs to load and sort them 
tiffs=dir(fullfile(tiffsPath, '*.tif'));
names = {tiffs.name};
[~,ndx] = natsort(names);
tiffs=tiffs(ndx);

% see what final size of image will be by loading a sample 
im = imread(fullfile(tiffs(1).folder,tiffs(1).name));
im = imresize(im, ResizeRatio, 'nearest');
[R,C] = size(im);
dual = false;
if R ~= C 
    R=min(size(im));
    C=R;
    dual = true;
end

if numFrames < 0
    n = length(tiffs);
else
    n = numFrames;
end

% get timestamp from first tiff file to find the number of initial frame drops and set up our drop detector
text = fileread(fullfile(tiffs(1).folder,tiffs(1).name));
k = strfind(text,'Time_From_Start = ');
timestamp = text(k+18:k+18+12);
prev_timestamp = str2double(timestamp(1:2))*60*60+str2double(timestamp(4:5))*60+str2double(timestamp(7:8))+str2double(timestamp(9:13));
offset = round(prev_timestamp/(1/fs));


% preallocate movie variable space... num intended frames = sum(frameDiffs)+offset+1
if numFrames < 0
    mov = nan(R,C,n+offset+1,'single');
else
    mov = nan(R,C,n,'single');
end


it = 1+offset;

if dual 
    % iterate through each file loading and performing relevant alterations
    for i =1:n
        text = fileread(fullfile(tiffs(i).folder,tiffs(i).name));
        k = strfind(text,'Time_From_Start = ');
        timestamp = text(k+18:k+18+12);
        timestamp = str2double(timestamp(1:2))*60*60+str2double(timestamp(4:5))*60+str2double(timestamp(7:8))+str2double(timestamp(9:13));

        im = imread(fullfile(tiffs(i).folder,tiffs(i).name));
        im = imresize(im, ResizeRatio, 'nearest');
        
        it = it + round((timestamp-prev_timestamp)/(1/fs));

        % based on i, determine which side of image to take
        if length(Order) == 3
            if rem(it, length(Order)) == 1
                image2take = Order(1);
            elseif rem(it, length(Order)) == 2
                image2take = Order(2);
            elseif rem(it, length(Order)) == 0
                image2take = Order(3);
            end
        elseif length(Order) == 4
            if rem(it, length(Order)) == 1
                image2take = Order(1);
            elseif rem(it, length(Order)) == 2
                image2take = Order(2);
            elseif rem(it, length(Order)) == 3
                image2take = Order(3);
            elseif rem(it, length(Order)) == 0
                image2take = Order(4);
            end
        elseif length(Order) == 2
            if rem(it, length(Order)) == 1
                image2take = Order(1);
            elseif rem(it, length(Order)) == 0
                image2take = Order(2);
            end
        end

        mov(:,:,it) = im(:,((image2take-1)*C+1):(image2take*C));
            
        prev_timestamp = timestamp;
        if mod(i,10000) == 0
            disp(['Read ' num2str(i) ' frames']);
        end
        
    end
else
    % iterate through each file loading and performing relevant alterations
    for i =1:n
        text = fileread(fullfile(tiffs(i).folder,tiffs(i).name));
        k = strfind(text,'Time_From_Start = ');
        timestamp = text(k+18:k+18+12);
        timestamp = str2double(timestamp(1:2))*60*60+str2double(timestamp(4:5))*60+str2double(timestamp(7:8))+str2double(timestamp(9:13));
        im = imread(fullfile(tiffs(i).folder,tiffs(i).name));
        im = imresize(im, ResizeRatio);
        it = it + round((timestamp-prev_timestamp)/(1/fs));
        mov(:,:,it) = im;
            
            prev_timestamp = timestamp;
        
    end
end

%close(wb);

end