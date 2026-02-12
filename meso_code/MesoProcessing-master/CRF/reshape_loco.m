function new_loco = reshape_loco(loco, t)


new = loco';

new_loco = nan(size(new, 1), length(t)-1, size(new, 2));

for i = 1:size(new, 2)
    test = new(:, i);
    test1 = repmat(test, 1, length(t)-1);
    new_loco(:,:,i) = test1;
end



end



