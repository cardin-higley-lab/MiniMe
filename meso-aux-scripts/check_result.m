main_folder = '/vast/palmer/scratch/higley/hd362/HD_Mouse_Training';

matfile = dir(fullfile(main_folder, '*Mouse2*', '*', 'Empty', 'meanvalue.mat'));

for i = 1:length(matfile)
    load(fullfile(matfile(i).folder, matfile(i).name))
    figure;
    plot(meanvalue)
    title(matfile(i).folder(end-22:end-18))
end
