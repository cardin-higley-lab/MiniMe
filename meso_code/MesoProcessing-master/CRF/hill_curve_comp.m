


%% load in loco data 

%WTs first
wt_loco_sessions = {'D:\meso_demo\CDKL5_CRF\ro130-339\339_crf_11122022\curve_data_loco_iter.mat',...
                'D:\meso_demo\CDKL5_CRF\ro139_365\365_crf_11112022\curve_data_loco_iter.mat',...
                'D:\meso_demo\CDKL5_CRF\ro139_365\365_crf_11142022\curve_data_loco_iter.mat',...
                'D:\meso_demo\CDKL5_CRF\ro139_365\365_crf_11152022\curve_data_loco_iter.mat',...
                'D:\meso_demo\CDKL5_CRF\ro139_365\365_crf_11162022\curve_data_loco_iter.mat',...
                'D:\meso_demo\CDKL5_CRF\ro159-406\406_01272023_crf\curve_data_loco_iter.mat'};
wt_loco_data = struct();
for i = 1:length(wt_loco_sessions)
    load(cell2mat(wt_loco_sessions(i)));
    wt_loco_data.c(i) = c;
    wt_loco_data.c50(i) = c50;
    wt_loco_data.rmax(i) = rmax;
    wt_loco_data.fitted_curve{i} = rmax./(1+((c50./(1:100)).^c));
end


% mutants next
mut_loco_sessions = {'D:\meso_demo\CDKL5_CRF\ro139_363\363_crf_11122022\curve_data_loco_iter.mat',...
                'D:\meso_demo\CDKL5_CRF\ro139_363\363_crf_11152022\curve_data_loco_iter.mat',...
                'D:\meso_demo\CDKL5_CRF\ro139_366\366_crf_11142022\curve_data_loco_iter.mat',...
                'D:\meso_demo\CDKL5_CRF\ro159-405\405_01262023_crf\curve_data_loco_iter.mat'};
mut_loco_data = struct();
for i = 1:length(mut_loco_sessions)
    load(cell2mat(mut_loco_sessions(i)));
    mut_loco_data.c(i) = c;
    mut_loco_data.c50(i) = c50;
    mut_loco_data.rmax(i) = rmax;
    mut_loco_data.fitted_curve{i} = rmax./(1+((c50./(1:100)).^c));
end


temp = nan(max(length(wt_loco_sessions), length(mut_loco_sessions)), 2);
temp(1:length(wt_loco_sessions), 1) = wt_loco_data.c50';
temp(1:length(mut_loco_sessions), 2) = mut_loco_data.c50';

figure(); boxchart(temp)



%% plot loco CRFs for 

wt_loco_avg = mean(cat(1, wt_loco_data.fitted_curve{:}), 1);
wt_loco_se = mean(cat(1, wt_loco_data.fitted_curve{:}), 1)./sqrt(length(wt_loco_sessions));
mut_loco_avg = mean(cat(1, mut_loco_data.fitted_curve{:}), 1);
mut_loco_se = mean(cat(1, mut_loco_data.fitted_curve{:}), 1)./sqrt(length(mut_loco_sessions));


hyper_curve = figure;hold on;

for i = 1:length(wt_loco_sessions)
    plot(wt_loco_data.fitted_curve{i},'k');
    xlim([0 100]); ylim([0 1.2])
end

for i = 1:length(mut_loco_sessions)
    plot(mut_loco_data.fitted_curve{i},'r');
    xlim([0 100]); ylim([0 1.2])
end


h1 = plot(1:100,wt_loco_avg, 'k.'); %plot actual dff values as dots 
shadedErrorBar(1:100, wt_loco_avg, wt_loco_se, 'lineprops', {'k', 'markerfacecolor','k'}); 
h2 = plot(1:100,mut_loco_avg, 'r.'); %plot actual dff values as dots 
shadedErrorBar(1:100, mut_loco_avg, mut_loco_se, 'lineprops', {'r', 'markerfacecolor','r'}); 

legend([h1(1) h2(1)], {'CDKL5 (+/Y)', 'CDKL5 (-/Y)'})





%% load in sit data 

%WTs first
wt_sit_sessions = {'D:\meso_demo\CDKL5_CRF\ro130-339\339_crf_11122022\curve_data_sit_iter.mat',...
                'D:\meso_demo\CDKL5_CRF\ro139_365\365_crf_11112022\curve_data_sit_iter.mat',...
                'D:\meso_demo\CDKL5_CRF\ro139_365\365_crf_11142022\curve_data_sit_iter.mat',...
                'D:\meso_demo\CDKL5_CRF\ro139_365\365_crf_11152022\curve_data_sit_iter.mat',...
                'D:\meso_demo\CDKL5_CRF\ro139_365\365_crf_11162022\curve_data_sit_iter.mat',...
                'D:\meso_demo\CDKL5_CRF\ro159-406\406_01272023_crf\curve_data_sit_iter.mat'};
wt_sit_data = struct();
for i = 1:length(wt_sit_sessions)
    load(cell2mat(wt_sit_sessions(i)));
    wt_sit_data.c(i) = c;
    wt_sit_data.c50(i) = c50;
    wt_sit_data.rmax(i) = rmax;
    wt_sit_data.fitted_curve{i} = rmax./(1+((c50./(1:100)).^c));
end


% mutants next
mut_sit_sessions = {'D:\meso_demo\CDKL5_CRF\ro139_363\363_crf_11122022\curve_data_sit_iter.mat',...
                'D:\meso_demo\CDKL5_CRF\ro139_363\363_crf_11152022\curve_data_sit_iter.mat',...
                'D:\meso_demo\CDKL5_CRF\ro139_366\366_crf_11142022\curve_data_sit_iter.mat',...
                'D:\meso_demo\CDKL5_CRF\ro159-405\405_01262023_crf\curve_data_sit_iter.mat'};
mut_sit_data = struct();
for i = 1:length(mut_sit_sessions)
    load(cell2mat(mut_sit_sessions(i)));
    mut_sit_data.c(i) = c;
    mut_sit_data.c50(i) = c50;
    mut_sit_data.rmax(i) = rmax;
    mut_sit_data.fitted_curve{i} = rmax./(1+((c50./(1:100)).^c));
end


temp = nan(max(length(wt_sit_sessions), length(mut_sit_sessions)), 2);
temp(1:length(wt_sit_sessions), 1) = wt_sit_data.c50';
temp(1:length(mut_sit_sessions), 2) = mut_sit_data.c50';

figure(); boxchart(temp)



%% plot sit CRFs for 

wt_sit_avg = mean(cat(1, wt_sit_data.fitted_curve{:}), 1);
wt_sit_se = mean(cat(1, wt_sit_data.fitted_curve{:}), 1)./sqrt(length(wt_sit_sessions));
mut_sit_avg = mean(cat(1, mut_sit_data.fitted_curve{:}), 1);
mut_sit_se = mean(cat(1, mut_sit_data.fitted_curve{:}), 1)./sqrt(length(mut_sit_sessions));


hyper_curve = figure;hold on;

for i = 1:length(wt_sit_sessions)
    plot(wt_sit_data.fitted_curve{i},'k');
    xlim([0 100]); ylim([0 1.2])
end

for i = 1:length(mut_sit_sessions)
    plot(mut_sit_data.fitted_curve{i},'r');
    xlim([0 100]); ylim([0 1.2])
end


h1 = plot(1:100,wt_sit_avg, 'k.'); %plot actual dff values as dots 
shadedErrorBar(1:100, wt_sit_avg, wt_sit_se, 'lineprops', {'k', 'markerfacecolor','k'}); 
h2 = plot(1:100,mut_sit_avg, 'r.'); %plot actual dff values as dots 
shadedErrorBar(1:100, mut_sit_avg, mut_sit_se, 'lineprops', {'r', 'markerfacecolor','r'}); 

legend([h1(1) h2(1)], {'CDKL5 (+/Y)', 'CDKL5 (-/Y)'})





































%% load in loco data 

%WTs first
wt_loco_sessions = {'D:\meso_demo\CRF_output\ro139_365_curve_data_loco.mat',...
                'D:\meso_demo\CRF_output\ro159-406_curve_data_loco.mat'};
wt_loco_data = struct();
for i = 1:length(wt_loco_sessions)
    load(cell2mat(wt_loco_sessions(i)));
    wt_loco_data.c(i) = c;
    wt_loco_data.c50(i) = c50;
    wt_loco_data.rmax(i) = rmax;
    wt_loco_data.fitted_curve{i} = fitted_curve;
end

% mutants next
mut_loco_sessions = {%'D:\meso_demo\CRF_output\ro139_363_curve_data_loco.mat',...
                'D:\meso_demo\CRF_output\ro139_366_curve_data_loco.mat',...
                'D:\meso_demo\CRF_output\ro159-405_curve_data_loco.mat'};
mut_loco_data = struct();
for i = 1:length(mut_loco_sessions)
    load(cell2mat(mut_loco_sessions(i)));
    mut_loco_data.c(i) = c;
    mut_loco_data.c50(i) = c50;
    mut_loco_data.rmax(i) = rmax;
    mut_loco_data.fitted_curve{i} = fitted_curve;
end


temp = nan(max(length(wt_loco_sessions), length(mut_loco_sessions)), 2);
temp(1:length(wt_loco_sessions), 1) = wt_loco_data.c50';
temp(1:length(mut_loco_sessions), 2) = mut_loco_data.c50';

figure(); boxchart(temp)
