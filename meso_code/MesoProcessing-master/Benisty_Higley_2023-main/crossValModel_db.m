function [Vm_tr, Vm, cBeta, cR, subIdx, cRidge, cLabels] =  crossValModel_db(fullR, Vc, cLabels, regIdx, regLabels, folds)
% function to compute cross-validated R^2
% made one change from original - circular permutation for
% folds.

cIdx = ismember(regIdx, find(ismember(regLabels,cLabels))); %get index for task regressors
cLabels = regLabels(sort(find(ismember(regLabels,cLabels)))); %make sure motorLabels is in the right order

%create new regressor index that matches motor labels
subIdx = regIdx;
subIdx = subIdx(cIdx);
temp = unique(subIdx);
for x = 1 : length(temp)
    subIdx(subIdx == temp(x)) = x;
end
cR = fullR(:,cIdx);

Vm = zeros(size(Vc),'single'); %pre-allocate motor-reconstructed V
s = RandStream('mt19937ar','Seed',1);% for reproducibility
foldCnt = floor(size(Vc,1) / folds);
randIdx = randi(s,foldCnt); %generate randum number index
cBeta = cell(1,folds);

for iFolds = 1:folds
    dataIdx = true(1,size(Vc,1));
    dataIdx(randIdx:randIdx+foldCnt) = false;
    
    if folds > 1
        dataIdx = circshift(dataIdx,randIdx+foldCnt*iFolds);
        %dataIdx(randIdx(((iFolds - 1)*foldCnt) + (1:foldCnt))) = false; %index for training data
        if iFolds == 1
            [cRidge, cBeta{iFolds}] = ridgeMML(Vc(dataIdx,:), cR(dataIdx,:), true); %get beta weights and ridge penalty for task only model
        else
            [~, cBeta{iFolds}] = ridgeMML(Vc(dataIdx,:), cR(dataIdx,:), true, cRidge); %get beta weights for task only model. ridge value should be the same as in the first run.
        end
        
        Vm(~dataIdx,:) = (cR(~dataIdx,:) * cBeta{iFolds}); %predict remaining data
        
%         if rem(iFolds,folds/5) == 0
%             fprintf(1, 'Current fold is %d out of %d\n', iFolds, folds);
%         end
    else
        [cRidge, cBeta{iFolds}] = ridgeMML(Vc, cR, true); %get beta weights for task-only model.
        Vm = (cR * cBeta{iFolds}); %predict remaining data
        disp('Ridgefold is <= 1, fit to complete dataset instead');
    end
end
% get train stats
[cRidge_tr, cBeta_tr] = ridgeMML(Vc, cR, true); %get beta w
Vm_tr = (cR * cBeta_tr); %predict remaining data
        
