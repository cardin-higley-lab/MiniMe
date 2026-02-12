function [res, cBeta, rejIdx, V] = train_test_ridge(Xz, Yz, Kfolds, inds2shuffle)

if exist('inds2shuffle', 'var') && ~isempty(inds2shuffle)
    sh_inds = randperm(size(Xz, 2));
    Xz(inds2shuffle, :) = Xz(inds2shuffle, sh_inds);
end
notnans = ~isnan(sum(Yz,2)) & ~isnan(sum(Xz, 1)');
Xz = Xz(:, notnans');
Yz = Yz(notnans, :);
Yz = zscore(Yz);
Xz = zscore(Xz')';
[~, fullQRR] = qr(bsxfun(@rdivide,Xz',sqrt(sum(Xz'.^2))),0); %orthogonalize normalized design matrix
rejIdx = false(1,size(Xz,1));

if sum(abs(diag(fullQRR)) > max(size(Xz')) * eps(fullQRR(1))) < size(Xz',2) %check if design matrix is full rank
    temp1 = ~(abs(diag(fullQRR)) > max(size(Xz')) * eps(fullQRR(1)));
    fprintf('Design matrix is rank-defficient. Removing %d/%d additional regressors.\n', sum(temp1), sum(~rejIdx));
    rejIdx(~rejIdx) = temp1; %reject regressors that cause rank-defficint matrix
end

Xz = Xz( ~rejIdx, :);

[V_tr,V, cBeta] = crossValModel_db(Xz', Yz, 1, ones(size(Xz, 1),1), 1, Kfolds);
for i = 1:size(Yz, 2)
    R_tr(i) = getRsq(double(V_tr(:,i)), Yz(:,i));
    R_te(i) = getRsq(double(V(:,i)), Yz(:,i));
end 

res.R_tr = max(R_tr, 0);
res.R_te = max(R_te, 0);

res.C_tr = diag(corr(V_tr, Yz))';
res.C_te = diag(corr(V, Yz))';