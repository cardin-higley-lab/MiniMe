function meso_clip(full_vid, filename, index, smoothed)


norm_clip = full_vid(:, index);

switch smoothed
    case 'smoothed_yes'
        norm_clip = smooth3(reshape(norm_clip, 256, 256, []));
    case 'smoothed_no'
        norm_clip = reshape(norm_clip, 256, 256, []);
end

norm_clip(:, 119:138, :) = nan;

for i = 1:size(norm_clip,3)
    h=figure('Position', [100 100 900 700]);
    toShow = norm_clip(:,:,i);
    h2 = imagesc(toShow);
    set(h2,'alphadata',~isnan(toShow));
    axis off
    low_lim = quantile(norm_clip(:), 0.01);
    high_lim = quantile(norm_clip(:), 0.99);
    clim([low_lim high_lim]);
    drawnow;
    frame = getframe(h);
    im = frame2im(frame);
    [imind,cm] = rgb2ind(im,256);
    if i == 1
          imwrite(imind,cm,filename, 'gif', 'Loopcount',inf,'DelayTime',.1); %
    else
          imwrite(imind,cm,filename, 'gif','WriteMode','append','DelayTime',.1); %
    end
    close(h)
end




end