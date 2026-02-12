function [tform] = register2channels(moving,fixed,R,C,outputPath)
%register two color images acquired with two cameras
minval=min(fixed,[],2);
maxval = max(fixed,[],2);
minval2=min(moving,[],2);
maxval2 = max(moving,[],2);


fixed_max=imadjust(mat2gray(reshape(maxval,R,C)));
moving_max=imadjust(mat2gray(reshape(maxval2,R,C)));
fixed_diffframe1 = reshape(maxval-minval, R, C);
moving_diffframe2 = reshape(maxval2-minval2, R, C);

%register images automatically first
try
    [MOVINGREG] = registerImages_1(moving_diffframe2,fixed_diffframe1);
    tform=MOVINGREG.Transformation;
    fixedRefObj=MOVINGREG.SpatialRefObj;
    movingRefObj=imref2d(size(moving_diffframe2));
    moving_reg = imwarp(moving_max, movingRefObj, tform, 'OutputView', fixedRefObj, 'Fillvalues',0);
    str='N';
    %while user is not satisfied with registration, keep asking if user wants
    %to try a different method for registration
    while strcmp(str,'N')
        %plot existing registration images
        h3=figure;subplot(2,2,1),imshow(fixed_max);title('Before Fixed');
        subplot(2,2,2),imshow(moving_max);title('Before Moving');
        subplot(2,2,3),imshow(fixed_max);title('After Fixed');
        subplot(2,2,4),imshow(moving_reg);title('AfterMoving');
        disp(num2str(tform.T));
        outFile=fullfile(outputPath,'GreenBlueRegistration.fig');
        saveas(h3,outFile,'fig');
        
        %create a flickering tif file to see before and after registration
        h1=fixed_max;
        outFile=fullfile(outputPath,'Before_reg.tif');
        imwrite(im2uint16(h1),outFile ,'tif');
        h1 = moving_max;
        imwrite(im2uint16(h1), outFile,'WriteMode', 'append');
        
        h2=fixed_max;
        outFile=fullfile(outputPath,'After_reg.tif');
        imwrite(im2uint16(h2),outFile ,'tif');
        h2 =moving_reg;
        imwrite(im2uint16(h2), outFile,'WriteMode', 'append');
        
        
        % check the output of registration amd allow interactive use if not good
        prompt = 'Is the registration satisfactory, if not we will do it interactively? Y/N [Y]: ';
        str = input(prompt,'s');
        if isempty(str)|| strcmp(str,'Y')
            disp('Continuing on without additional modification to registration')
            break;
        end
        clear tform MOVINGREG;evalin( 'base', 'clear movingReg' )
        prompt = 'Automated interactive(A) or control point registration(C)? A/C [C]: ';
        string1= input(prompt,'s');
        if isempty(string1)|| strcmp(string1,'C')
            [movingPoints, fixedPoints]=cpselect(moving_max,fixed_max,'Wait',true);%do control point selection with a gui
            tform = fitgeotrans(movingPoints,fixedPoints,'NonReflectiveSimilarity');%no option of rigid transformation so use NonReflectiveSimilarity instead
            moving_reg = imwarp(moving_max,tform,'OutputView',imref2d(size(fixed_diffframe1)));
        elseif strcmp(string1,'A')
            disp('This will open the matlab GUI for registration estimation. Start with multimodal intensity registration with followin settings: Normalize: Yes; Transformation: rigid')
            disp('Use default parameters but can be helpful to reduce radius to increase registration quality');
            disp('When done, click export to export to workspace');
            registrationEstimator(moving_diffframe2,fixed_diffframe1);
            prompt = 'Press any key when you are done: ';
            string2 = input(prompt,'s');
            MOVINGREG = evalin('base','movingReg');
            tform=MOVINGREG.Transformation;
            fixedRefObj=MOVINGREG.SpatialRefObj;
            movingRefObj=imref2d(size(moving_diffframe2));
            moving_reg = imwarp(moving_max, movingRefObj, tform, 'OutputView', fixedRefObj, 'Fillvalues',0);
        end
    end
catch % if registration fails to converge do manually
       disp('Registration Error. Do manually using control points'); 
        [movingPoints, fixedPoints]=cpselect(moving_max,fixed_max,'Wait',true);%do control point selection with a gui
        tform = fitgeotrans(movingPoints,fixedPoints,'NonReflectiveSimilarity');%no option of rigid transformation so use NonReflectiveSimilarity instead
        moving_reg = imwarp(moving_max,tform,'OutputView',imref2d(size(fixed_diffframe1)));
      

        %plot registration images
        h3=figure;subplot(2,2,1),imshow(fixed_max);title('Before Fixed');
        subplot(2,2,2),imshow(moving_max);title('Before Moving');
        subplot(2,2,3),imshow(fixed_max);title('After Fixed');
        subplot(2,2,4),imshow(moving_reg);title('AfterMoving');
        disp(num2str(tform.T));
        outFile=fullfile(outputPath,'ChannelRegistration.fig');
        saveas(h3,outFile,'fig');
        
        %create a flickering gif file to see before and after registration
        filename = fullfile(outputPath, 'Before_reg.gif');
        before_im = cat(3, fixed_max, moving_max);
        for i = 1:2
            h=figure();
            toShow = before_im(:,:,i);
            h2 = imagesc(toShow);
            set(h2,'alphadata',~isnan(toShow));
            axis off
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


        %create a flickering gif file to see before and after registration
        filename = fullfile(outputPath, 'After_reg.gif');
        after_im = cat(3, fixed_max, moving_reg);
        for i = 1:2
            h=figure();
            toShow = after_im(:,:,i);
            h2 = imagesc(toShow);
            set(h2,'alphadata',~isnan(toShow));
            axis off
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

end




