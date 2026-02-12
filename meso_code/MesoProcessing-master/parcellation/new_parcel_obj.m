%% let's try to divide up larger allen parcels

%% load some data to play w later
load('D:\meso_demo\CDKL5\ro096-216\216_06132022\Ca_traces_spt_patch9_Allen.mat')
allen_traces = parcels_time_trace;
numAllen = size(allen_traces, 1); 

%% load parcel obj and plot
load('parcel_Obj.mat');

% make and show mapping 
allen_map = nan(256,256);
allen_sizes = nan(size(parcel_obj.indicators, 3), 2);
for i = [1:22 27:56]
    allen_sizes(i, :) = [i, sum(parcel_obj.indicators(:,:,i)>0, 'all')];
    allen_map(parcel_obj.indicators(:,:,i)>0)=i;
end
figure; imagesc(allen_map);

%% define parcels that need to be broken up and plot them 
long_allen = [17:22, 49:54];
for i = 1:2:length(long_allen)
    parcel_map = parcel_obj.indicators(:,:,long_allen(i)); %parcelNum
    figure; imagesc(parcel_map);
end


%% make new indicators
%numParcels = size(parcel_obj.indicators, 3);
%indicators = nan(256, 256, numParcels+long_allen);

numNewInds = numParcels+long_allen;

indicators = [];
names = {};

for i = 1:size(parcel_obj.indicators, 3)

    if ~any(i == long_allen) 
        indicators = cat(3, indicators, parcel_obj.indicators(:,:,i));
        names(end+1) = parcel_obj.names(i);
    
    else  
        parcel_size = allen_sizes(i, 2); 
        new_size = round(parcel_size/2);
    
        parcel_map = parcel_obj.indicators(:,:,i);
        numPixels = nan(size(parcel_map, 2), 1);
        for row = 1:length(numPixels)
            numPixels(row) = sum(parcel_map(1:row, :), 'all');
        end
    
        midway_row = find(numPixels>new_size, 1, 'first');
    
        new_parcel_map1 = [parcel_map(1:midway_row, :); zeros(size(parcel_map(midway_row+1:end, :)))];
        new_parcel_map2 = [zeros(size(parcel_map(1:midway_row, :))); parcel_map(midway_row+1:end, :)];
    
        indicators = cat(3, indicators, new_parcel_map1, new_parcel_map2);

        names(end+1) = {[char(parcel_obj.names(i)), '-posterior']};
        names(end+1) = {[char(parcel_obj.names(i)), '-anterior']};

    end


end

%% plot new map


new_allen_map = nan(256,256);
for i = [1:28 33:68]
    new_allen_map(indicators(:,:,i)>0)=i;
end
figure; imagesc(new_allen_map);



%% get new allen centers

numNewAllen = size(indicators, 3);

new_allen_centers = nan(numNewAllen, 2);
for i = 1:numNewAllen
    [r, c] = find(new_allen_map == i); 
    new_allen_centers(i, :) = [mean(r), mean(c)];
end



%% load in some lssc data 

load('D:\meso_demo\CDKL5\ro096-216\216_06132022\Ca_traces_spt_patch9_LSSC.mat')
lssc_traces = parcels_time_trace; 
numLSSC = size(lssc_traces, 1); % get the number of lssc parcels
lssc_map = parmap;


%get centers
lssc_centers = nan(numLSSC, 2);
lssc_sizes = nan(numLSSC, 2);
for i = 1:numLSSC
    [r, c] = find(lssc_map == i); 
    lssc_sizes(i, :) = [i, length(find(lssc_map == i))];
    lssc_centers(i, :) = [mean(r), mean(c)];
end

% remove small lssc parcels 
% ordered_lssc = sortrows(lssc_sizes, 2);
% smallest_lssc = find(lssc_sizes(:,2) < 200);
smallest_lssc = [44, 46, 88, 92];
lssc_centers(smallest_lssc,:) = NaN;


% adjust lssc map too
for i = length(smallest_lssc)
    lssc_map(lssc_map==smallest_lssc(i)) = NaN;
end



%% make matches... 
  

distances = zeros(numLSSC, numNewAllen);
for lssc = 1:numLSSC
    for allen = 1:numNewAllen
        lssc_coord = lssc_centers(lssc, :);
        allen_coord = new_allen_centers(allen, :);
        distances(lssc,allen) = pdist([lssc_coord; allen_coord], 'euclidean');
    end 
end


% choose the closest lssc center for each allen center
% smallest_dists = min(distances);
% [row, col] = find(distances==smallest_dists);
% pairs = [row, col];



% more advanced way to make pairing after finding distances 

smallest_dists = nan(numNewAllen, 3, numNewAllen);
for i = 1:numNewAllen
    vector = [i:numNewAllen, 1:i-1];
    dist_vec = distances;
    for allen = vector
        [dist, lssc_num] = min(dist_vec(:,allen));
        if isnan(dist)
            smallest_dists(allen, :, i) = [allen, NaN, NaN];
        else
            smallest_dists(allen, :, i) = [allen, lssc_num, dist];
            dist_vec(lssc_num, :) = NaN;
        end
    end
end



cum_dists = nan(numNewAllen, 1);
for i = 1:numNewAllen
    cum_dists(i) = sum(smallest_dists(:,3,i), 'omitnan');
end

[best_solution, index] = min(cum_dists);

smallest_dists = smallest_dists(:,:, index);


%% Plot

% plot figs to show relationship 
temp = nan(256, 256);
for i = 1:length(smallest_dists)
    temp(lssc_map==smallest_dists(i, 2)) = i;
end



figure(); tiledlayout(1, 3);
nexttile(); imagesc(lssc_map+1); 
nexttile(); imagesc(new_allen_map); caxis([0 68])
nexttile(); imagesc(temp); caxis([0 68])
colormap([ 1 1 1; parula])




% % plot figs to show relationship 
% temp = nan(256, 256);
% for i = 1:length(pairs)
%     temp(lssc_map==pairs(i, 1)) = i;
% end
% 
% 
% 
% figure(); tiledlayout(1, 3);
% nexttile(); imagesc(lssc_map+1); 
% nexttile(); imagesc(new_allen_map); caxis([0 68])
% nexttile(); imagesc(temp); caxis([0 68])
% colormap([ 1 1 1; parula])

%% problems and future directions 

length(unique(pairs(:,1)))

% so above signals that not all allen parcels are represented by a unique
% lssc, but some lsscs are representative of more than one allen 

% to optimize, thought about looping through distances to find th shortest
% and removing the chosen parcel from the list of options so that there
% isnt a single lssc that can be chosen for more than one allen ... ideally
% would do this iteratively, starting with each consecutive allen parcel to
% see which matching best minimizes distances overall 



% hadas's take: ditch the allen entirely
% [1:26 PM] for non-correlation matrix analyses: average functional maps across animals, 
% pixelwise and can then overlay allen at the end if necessary to talk about regions
% [1:27 PM] for correlation matrices, she suggests using a grid system to identify specific 
% functional parcels, rather than the allen, which will give a better distribution and more coverage across the brain
% [1:28 PM] she did show me how to adjust the lssc to give fewer parcels so i can still do that and see 
% how it works with the code i have for converting lssc to allen if we're curious, but if you both agree with her take, then ill ditch it
% 3:37
% mike also mentioned that for the correlation matrices... can do corr mats using lssc and then collapse on the 
% row to get the corr-value for each lssc parcel... can then permute the ones that fall within v1 for ex


