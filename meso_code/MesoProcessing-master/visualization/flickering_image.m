function flickering_image(image1, image2, filename)

after_im = cat(3, image1, image2);
for i = 1:2
    h=figure();
    toShow = after_im(:,:,i);
    h2 = imagesc(toShow);
    set(h2,'alphadata',~isnan(toShow));
    clim([quantile(image1, 0.05, 'all') quantile(image1, 0.95, 'all')])
    axis off
    drawnow;
    frame = getframe(h); im = frame2im(frame);
    [imind,cm] = rgb2ind(im,256);
    if i == 1
        imwrite(imind,cm,filename, 'gif', 'Loopcount',inf,'DelayTime',.5); %
    else
        imwrite(imind,cm,filename, 'gif','WriteMode','append','DelayTime',.5); %
    end
    close(h)
end

end