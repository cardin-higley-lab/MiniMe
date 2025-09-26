function saveColumnwise(data,path,prefix) %% save
mkdir(path)
for i = 1:size(data,2)   
    column = squeeze(data(:,i,:));
    save(fullfile(path,[prefix 'Col' num2str(i) '.mat']),"column");
end
