function [dff_parcels, parcellationmap] = functional_parcellation_LSSC(regsig, cfg, brainMask)


brain_mask_left = brainMask;
brain_mask_right = brainMask;
brain_mask_left(:,129:end) = false; %masks out all pixels on the right
brain_mask_right(:,1:128) = false; %masks out all pixels on the left

allregionspix_left = find(brain_mask_left);
allregionspix_right = find(brain_mask_right);

[R,C] = size(brainMask);
dFoF_masked_left = regsig(allregionspix_left, :);
dFoF_masked_right = regsig(allregionspix_right, :);

cfg.preProcess=false;
cfg.N_TRIALS=10;        %trials = segments of the session (https://www.biorxiv.org/content/10.1101/2021.08.15.456390v1.full.pdf)
cfg.n_clust = 5000 ;
cfg.makePlots = true;
cfg.thrcluster=0.99;
cfg.NROWS = R;
cfg.NCOLS = C;
cfg.isoverlap = false;

% first we process the left hemisphere
cfg.title_str = 'left';
leftfile = dir(fullfile(cfg.outputfilepath, 'left*.mat')); % looks to see if there's a *left file already
if isempty(leftfile) % if no left file (hasnt processed), then proceeds 
    runROI_meso_nlm(cfg, dFoF_masked_left, allregionspix_left, brain_mask_left);
    leftfile = dir(fullfile(cfg.outputfilepath, 'left*.mat'));
end
leftparcels = load(fullfile(cfg.outputfilepath,leftfile.name));


% next we process the right hemisphere
cfg.title_str = 'right';
runROI_meso_nlm(cfg, dFoF_masked_right, allregionspix_right, brain_mask_right);
rightfile = dir(fullfile(cfg.outputfilepath, 'right*.mat'));
if isempty(rightfile)
    runROI_meso_nlm(cfg, dFoF_masked_left, allregionspix_left, brain_mask_left);
    rightfile = dir(fullfile(cfg.outputfilepath, 'right*.mat'));
end
rightparcels = load(fullfile(cfg.outputfilepath,rightfile.name));

% then you do something else... 
ROI_list = [leftparcels.ROI_list  rightparcels.ROI_list];
dff_parcels = nan(length(ROI_list), size(regsig,2));
for k = 1:length(ROI_list)
    parcellationmap(ROI_list(k).pixel_list) = k;
    dff_parcels(k, :) = nanmean(regsig(ROI_list(k).pixel_list, :));
end
