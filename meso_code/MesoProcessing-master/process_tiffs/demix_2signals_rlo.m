function [signal1, signal2, skippedframes, skippedchannels] = demix_2signals_rlo(mov, paramsSignalsExtraction)

ratio = paramsSignalsExtraction.blueuvRatio;

% assigns first and second frame for rest of channel alignments 
firstframe=mov(:,:,1);
secondframe=mov(:,:,2);

% takes zoomed in frames and linearizes for correlations later
firstframe_linear= double(reshape(firstframe,[],1));  %double(reshape(firstframe(68:188,68:188),[],1));
secondframe_linear= double(reshape(secondframe,[],1)); %double(reshape(secondframe(68:188,68:188),[],1));

skippedcount=0;           skippedframes=NaN(100,1);         skippedchannels=NaN(100,1); 

ch=0;                     ind_bl=1;                         ind_uv=1;

[R,C,nframes] = size(mov);
signal1=zeros(R*C,nframes+100,'uint16');
signal2=zeros(R*C,nframes+100,'uint16');

for iframe=1:nframes 
   
    if mod(iframe,2000)==0
        disp(['Checking Dropped Frames, done ',num2str(iframe),' frames']);
    end
    
    if iframe<10 % if within the first 10 frames... set blue and uv to 1st and 2nd frame
        blframe_linear=firstframe_linear;
        uvframe_linear=secondframe_linear; 
    
    elseif ind_uv-1 > 0 && ind_bl-1 > 0 %for all frames after the first two.. 
        x = reshape(signal1(:,ind_bl-1),R,C);
        blframe_linear = double(reshape(x,[],1));
        
        x = reshape(signal2(:,ind_uv-1),R,C);
        uvframe_linear = double(reshape(x,[],1));
    end
    
    currentframe = mov(:,:,iframe); %get current frame from mov
    currentframe_linear = double(reshape(currentframe,[],1));

    if corr(currentframe_linear, uvframe_linear) > corr(currentframe_linear, blframe_linear)
        ch_current=0; %if current frame is closer to uv... 
    else
        ch_current=1;
    end
    
    if ch_current==ch % frame drop
        skippedcount=skippedcount+1;
        display(['Dropped frame at ',num2str(iframe)]);
        skippedchannels(skippedcount)=~ch;
        skippedframes(skippedcount)=iframe;
        switch ch_current
            case 0 % skipped a blue
                signal1(:,ind_bl)=signal1(:,max(ind_bl-1,1));
                signal2(:,ind_uv)=currentframe(:);
            case 1 % skipped a uv
                signal1(:,ind_bl)=currentframe(:);
                signal2(:,ind_uv)=signal2(:,max(ind_uv-1,2));
        end
        ind_bl=ind_bl+1; ind_uv=ind_uv+1;
        
    else  % if alternating properly
        switch ch_current
            case 0
                signal2(:,ind_uv)=currentframe(:);
                ind_uv=ind_uv+1;
            case 1
                signal1(:,ind_bl)=currentframe(:);
                ind_bl=ind_bl+1;
        end
        ch=~ch;
    end
    
    
end

skippedframes=skippedframes(1:skippedcount);
skippedchannels=skippedchannels(1:skippedcount); % record skipped frames ind and channels
signal1=signal1(:,1:ind_bl-1);
signal2=signal2(:,1:ind_uv-1);

skippedframes_asString=[];
for i = 1:length(skippedframes)
   skippedframes_asString = [skippedframes_asString, num2str(skippedframes(i)), ', '];
end


if ~isempty(skippedframes) == 0
    disp('no skipped frames')
else 
    disp(['skipped frames: ' skippedframes_asString]);
end

signal1 = single(signal1);
signal2 = single(signal2);

if ratio > 1
    signal1 = interpolate_pixels(signal1, ratio);
    signal2 = interpolate_pixels(signal2, ratio);
end

