%% compare lssc to allen


%% load lssc data, map, and calculate centers  

load('D:\meso_demo\CDKL5\ro096-216\216_06132022\Ca_traces_spt_patch9_LSSC.mat')
lssc_traces = parcels_time_trace; 
numLSSC = size(lssc_traces, 1); % get the number of lssc parcels
lssc_map = parmap;
clear parcels_time_trace parcels_time_trace_GSR parmap

% show mapping 
figure; imagesc(lssc_map)

%get centers
lssc_centers = nan(numLSSC, 2);
lssc_sizes = nan(numLSSC, 2);
for i = 1:numLSSC
    [r, c] = find(lssc_map == i); 
    lssc_sizes(i, :) = [i, length(find(lssc_map == i))];
    lssc_centers(i, :) = [mean(r), mean(c)];
end

%%  QUESTION    : should we manually remove certain parcels? ie 44 here is the thin blue around the top... 


% % manually remove small/weird lssc parcels 
% for i = 1:numLSSC
%     figure(); imagesc(lssc_map==i)
%     title(sprintf('parcel num: %d', i))
% end
% 
% parcels2remove = [92 88 46 44];
% 
% lssc_centers(parcels2remove,:) = NaN;
% 
% 
% % are these just the smallest parcels...?
% ordered_lssc = sortrows(lssc_sizes, 2);
% 
% smallest_lssc = ordered_lssc(1:8, 1); 
% for i = 1:length(smallest_lssc)
%     figure(); imagesc(lssc_map==smallest_lssc(i))
%     title(sprintf('parcel num: %d', smallest_lssc(i)))
% end


% programatically remove smallest parcels (only a few)
ordered_lssc = sortrows(lssc_sizes, 2);
smallest_lssc = find(lssc_sizes(:,2) < 200);
lssc_centers(smallest_lssc,:) = NaN;

% plot smallest parcels
for i = 1:length(smallest_lssc)
    figure(); imagesc(lssc_map==2)
end


% take 56 largest parcels 
smallest_lssc = ordered_lssc(1:end-55, 1);
lssc_centers(smallest_lssc,:) = NaN;


%% load parcel object and make allen parmap

load('D:\meso_demo\CDKL5\ro096-216\216_06132022\Ca_traces_spt_patch9_Allen.mat')
allen_traces = parcels_time_trace;
numAllen = size(allen_traces, 1); % get the number of lssc parcels 
load('parcel_Obj.mat');
clear parcels_time_trace parcels_time_trace_GSR

% make and show mapping 
allen_map = nan(256,256);
for i = [1:22 27:56]
    allen_map(parcel_obj.indicators(:,:,i)>0)=i;
end
figure; imagesc(allen_map);

% get centers
allen_centers = nan(numAllen, 2);
for i = 1:numAllen
    [r, c] = find(allen_map == i); 
    allen_centers(i, :) = [mean(r), mean(c)];
end


%% make matches... 
  

distances = zeros(numLSSC, numAllen);
for lssc = 1:numLSSC
    for allen = 1:numAllen
        lssc_coord = lssc_centers(lssc, :);
        allen_coord = allen_centers(allen, :);
        distances(lssc,allen) = pdist([lssc_coord; allen_coord], 'euclidean');
    end 
end


% choose the closest lssc center for each allen center
smallest_dists = min(distances);
[row, col] = find(distances==smallest_dists);


% row values are lssc parcel # that best correspond to allen parcel
pairs = [row, col];


% plot figs to show relationship 
figure();
temp = nan(256, 256);
for i = 1:length(pairs)
    temp(lssc_map==pairs(i, 1)) = i;
end
imagesc(temp)



figure; imagesc(allen_map);




figure(); tiledlayout(1, 3);
nexttile(); imagesc(lssc_map);
nexttile(); imagesc(allen_map);
nexttile(); imagesc(temp)


%% what about just mapping lssc to allen based on whether center of mass is within the allen parcel... therefore can map >1 lssc to each allen parcel


% 1. get lssc centers... still going to remove the smallest ones...
ordered_lssc = sortrows(lssc_sizes, 2);
smallest_lssc = find(lssc_sizes(:,2) < 200);
lssc_centers(smallest_lssc,:) = NaN;

% 2. map these centers to allen mapping... 

mappings = nan(numLSSC, 1);
for i = 1:numLSSC
    lssc_center = round(lssc_centers(i,:));
    if ~isnan(lssc_center(1))
        allen_parcel = allen_map(lssc_center(1), lssc_center(2));
        mappings(i) = allen_parcel;
    else
        mappings(i) = NaN;
    end
end



% plot figs to show relationship 
figure();
temp = nan(256, 256);
for i = 1:length(mappings)
    temp(lssc_map==i) = mappings(i);
end
imagesc(temp)


numUnique_lssc = length(unique(mappings)); % only maps to 42 allen (ie 14 allen did not contain a center of mass of an lssc parcel)




figure; imagesc(allen_map);





figure(); tiledlayout(1, 3);
nexttile(); imagesc(lssc_map);
nexttile(); imagesc(allen_map);
nexttile(); imagesc(temp)




% %% let's try to divide up larger allen parcels
% 
% 
% load('D:\meso_demo\CDKL5\ro096-216\216_06132022\Ca_traces_spt_patch9_Allen.mat')
% allen_traces = parcels_time_trace;
% numAllen = size(allen_traces, 1); % get the number of lssc parcels 
% load('parcel_Obj.mat');
% 
% 
% % make and show mapping 
% allen_map = nan(256,256);
% allen_sizes = nan(size(parcel_obj.indicators, 3), 2);
% for i = [1:22 27:56]
%     allen_sizes(i, :) = [i, sum(parcel_obj.indicators(:,:,i)>0, 'all')];
%     allen_map(parcel_obj.indicators(:,:,i)>0)=i;
% end
% figure; imagesc(allen_map);
% 
% 
% % 
% % ordered_allen = sortrows(allen_sizes, 2, 'descend');
% % ordered_allen = ordered_allen(~isnan(ordered_allen(:,1)), :);
% % largest_allen = ordered_allen(1:12, :);
% % [parcelSize, parcelNum] = max(ordered_allen(:, 2));
% % 
% % for i = 1:2:length(largest_allen)
% %     parcel_map = parcel_obj.indicators(:,:,largest_allen(i)); %parcelNum
% %     figure; imagesc(parcel_map);
% % end
% 
% 
% long_allen = [17:22, 49:54];
% for i = 1:2:length(long_allen)
%     parcel_map = parcel_obj.indicators(:,:,long_allen(i)); %parcelNum
%     figure; imagesc(parcel_map);
% end
% 
% indicators = parcel_obj.indicators;
% for i = 1:length(long_allen)
%     %i=1;
%     parcel = long_allen(i);
%     parcel_size = allen_sizes(parcel, 2); 
%     new_size = round(parcel_size/2);
% 
%     parcel_map = parcel_obj.indicators(:,:,parcel);
%     numPixels = nan(size(parcel_map, 2), 1);
%     for row = 1:length(numPixels)
%         numPixels(row) = sum(parcel_map(1:row, :), 'all');
%     end
% 
%     midway_row = find(numPixels>new_size, 1, 'first');
% 
%     new_parcel_map1 = [parcel_map(1:midway_row, :); zeros(size(parcel_map(midway_row+1:end, :)))];
%     new_parcel_map2 = [zeros(size(parcel_map(1:midway_row, :))); parcel_map(midway_row+1:end, :)];
% 
%     %% 1. parcel_obj.combinedParcels
%     % generally ignoring this bc doesnt actually seem super useful due to the boudaries 
% 
% %     temp = parcel_obj.combinedParcels;
% %     temp(new_parcel_map2==1) = parcel + 1;
% %     figure; imagesc(temp);    
% 
%     %% 2. parcel_obj.indicators
% 
%     %temp = parcel_obj.indicators;
%     indicators = cat(3, indicators(:,:,1:parcel-1), new_parcel_map1, new_parcel_map2, indicators(:,:,(parcel+1):end));
% 
%     %% 3. parcel_obj.names
% 
% %     parcel_obj.names
% % 
% %     new_names = 
% 
% end
% 
% 
% new_allen_map = nan(256,256);
% for i = 1:size(indicators, 3)
%     new_allen_map(indicators(:,:,i)>0)=i;
% end
% figure; imagesc(new_allen_map);

