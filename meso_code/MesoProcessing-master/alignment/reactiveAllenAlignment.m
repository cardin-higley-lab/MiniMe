function reactiveAllenAlignment(movingImage)
    %%reactiveAllenAlignment A reactive GUI for aligning data to allen atlas.
    % NOTE: THIS WILL SET/OVERWRITE THE VARIABLE 'tform' IN YOUR WORKSPACE.
    % This function used callbacks to update the transformed image in real-time
    % with alignment point movement, giving the user instant feedback on the
    % alignment quality.

    % load allen atlas
    load('parcells_updated121519.mat','parcells_new');
    % mask out colliculi
    mask = true(size(parcells_new.indicators,3),1);
    mask(23:26) = false;
    % get borders and flip to preferred orientation
    %fixedImage = rot90(edge(sum(parcells_new.indicators(:,:,mask),3)),2);
    fixedImage = rot90(edge(sum(parcells_new.indicators(:,:,mask),3)),2);
    
    % create generic fixed and moving points
    %fixedPoints = [round(size(fixedImage,1)/4) round(size(fixedImage,2)/4); round(3*size(fixedImage,1)/4) round(size(fixedImage,2)/4); round(size(fixedImage,1)/4) round(3*size(fixedImage,2)/4)];
    
    fixedPoints = [70.8576   62.6474; 187.8377   62.6474; 129.0651  184.9211];
    movingPoints = [round(size(movingImage,2)/4) round(size(movingImage,1)/4); round(3*size(movingImage,2)/4) round(size(movingImage,1)/4); round(size(movingImage,2)/4) round(3*size(movingImage,1)/4)];
    
    % generate initial tform
    tform = fitgeotform2d(movingPoints,fixedPoints,"similarity");

    fixedPointHandles = {};
    movingPointHandles = {};
    
    colors = {'r','g','b'};
    f = figure(42);
    % put fixed and moving image in location that will be easy to access
    % during callback function
    f.UserData = {fixedImage,movingImage};
    tiledlayout(1,3,'TileSpacing','tight');
    % first pane
    nexttile;
    imagesc(fixedImage)
    axis off
    title('Fixed')
    
    % add alignment points for fixed image
    for i =1:3
        fixedPointHandles{end+1} = images.roi.Point(gca,'Position',[fixedPoints(i,:)],'Color',colors{i},'UserData',['f' num2str(i)]);
        % connect point moving to function updating tform and transformed
        % image
        addlistener(fixedPointHandles{end},'MovingROI',@alignmentPointMoveEvent);
        addlistener(fixedPointHandles{end},'ROIMoved',@alignmentPointMoveEvent);
    end
    
    % second pane
    nexttile;
    imagesc(movingImage)
    axis off;
    title('Moving')

    % add alignment points for moving image
    for i =1:3
        movingPointHandles{end+1} = images.roi.Point(gca,'Position',[movingPoints(i,:)],'Color',colors{i},'UserData',['m' num2str(i)]);
        % connect point moving to functoin updating tform and transformed
        % image
        addlistener(movingPointHandles{end},'MovingROI',@alignmentPointMoveEvent);
        addlistener(movingPointHandles{end},'ROIMoved',@alignmentPointMoveEvent);
    end
    
    % third pane
    nexttile;
    
    % transform moving image by tform
    transformed = imwarp(movingImage,tform,OutputView=imref2d(size(fixedImage)));
    % display fused image
    transformHandle = imagesc(imfuse(transformed,fixedImage));
    axis off
    % mark this image to be update for later (makes it more convenient to
    % grab this image during callback function)
    transformHandle.UserData = true;
    title('Transformed')
    linkdata on
    
    % while loop to prevent program advancement until the figure is closed
    while size(findobj(f))>0
       pause(0.1)
    end

end

 function alignmentPointMoveEvent(src,evt)
    tp = ancestor(src,'figure','toplevel');
    points = findobj('Type','images.roi.Point');
    transformedHandle = findobj('UserData',true,'Type','image');
    movingPoints = zeros(3,2);
    fixedPoints = zeros(3,2);
    for p =1:length(points)
        if strcmp(points(p).UserData(1),'m')
            movingPoints(str2num(points(p).UserData(2)),:) = points(p).Position;
        else
            fixedPoints(str2num(points(p).UserData(2)),:) = points(p).Position;
        end
    end
    tform = fitgeotform2d(movingPoints,fixedPoints,"similarity");
    transformed = imwarp(tp.UserData{2},tform,OutputView=imref2d(size(tp.UserData{1})));
    assignin('base', 'tform', tform );
    assignin('base', 'transformed', transformed);
    set(transformedHandle,'Cdata',imfuse(transformed,tp.UserData{1}));
end