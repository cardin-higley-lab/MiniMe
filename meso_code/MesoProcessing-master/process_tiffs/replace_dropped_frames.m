function [blue_raw, uv_raw] = replace_dropped_frames(blue_raw, uv_raw)

%% replace frame drops in blue channel
blue_nans = find(squeeze(sum(isnan(blue_raw),[1 2]))==size(blue_raw,1)*size(blue_raw,2));
for it = 1:length(blue_nans)
    blue_raw(:,:,blue_nans(it)) = blue_raw(:,:,blue_nans(it)-1);
end

skippedframes_asString=[];
for i = 1:length(blue_nans)
   skippedframes_asString = [skippedframes_asString, num2str(blue_nans(i)), ', '];
end

if ~isempty(blue_nans)
    disp(['number of blue frames dropped: ' num2str(length(blue_nans))]);
    disp(['skipped frames: ' skippedframes_asString]);
end



%% replace frame drops in uv channel
uv_nans = find(squeeze(sum(isnan(uv_raw),[1 2]))==size(uv_raw,1)*size(uv_raw,2));
for it = 1:length(uv_nans)
    uv_raw(:,:,uv_nans(it)) = uv_raw(:,:,uv_nans(it)-1);
end

skippedframes_asString=[];
for i = 1:length(uv_nans)
   skippedframes_asString = [skippedframes_asString, num2str(uv_nans(i)), ', '];
end

if ~isempty(uv_nans)
    disp(['number of uv frames dropped: ' num2str(length(uv_nans))]);
    disp(['skipped frames: ' skippedframes_asString]);
end





end 