
clear vars
project = 'CDKL5_CRF';
dir_path = ['D:\meso_demo\' project];
masterDir = dir(dir_path);

for animals = 1:length(masterDir)

    animal = masterDir(animals).name;
    if length(animal) < 3
        continue;
    end

    %% combine within animal data to increase trials/animal
    loco_data = [];
    sit_data = [];

    % get subfolder data to added to aggregated data       
    subfolders = dir(fullfile(dir_path, animal));
    for sessions = 1:length(subfolders)

        session = subfolders(sessions).name;
        if length(session) < 3
            continue;
        end

        folder = fullfile(dir_path, animal, session);
        disp(['working on: ' folder])

            %% hyperbolic curve fitting for crf data
            
            load(fullfile(folder, 'crf_data_lowface10p.mat'));

            loco_responses = nan(size(v1_responses));
            loco_responses(loco==1) = v1_responses(loco==1);
            loco_data = [loco_data loco_responses];

            sit_responses = nan(size(v1_responses));
            sit_responses(loco==0) = v1_responses(loco==0);
            sit_data = [sit_data sit_responses];

    end

            contrasts = [1 2 5 10 20 50 100];
            x=contrasts';

            y = loco_data;
            yBS = zeros(size(y,1),1); %new array to hold bootstrap data
            yAvg = mean(y,2, 'omitnan');  %generate average of all data to compare with bootstrap
            iter = 10000;   %number of times to sample with replacement for bootstrap

            % run the bootstrap to produce yBS
            for i = 1:length(contrasts) %note - assumes data are (contrasts,trials)
                temp=y(i,:);
                temp=temp(isnan(temp)==0); %only pull out non-Nan points
                val=0;
                for j=1:iter
                    val=val+(temp(randi([1 length(temp)])));
                end
                yBS(i)=val/iter;
            end

            % curve fit to Hill equation
            fitfun = fittype(@(a,b,c,x) a./(1+((b./x).^c)));  % a=Rmax, b=c50, c=slope
            lower = [0,0,0]; %lower bound for params
            upper = [5,5000,10];  %upper bound for params
            start = [0.5,5,1];  %initial values for params
            [fitted_curve,gof] = fit(x,yBS,fitfun, 'Lower', lower, 'Upper', upper, 'StartPoint', start);
            %[fitted_curve,gof] = fit(x,yBS,fitfun, 'StartPoint', start);
            coeffs = coeffvalues(fitted_curve);
            rmax = coeffs(1); disp(num2str(rmax))
            c50 = coeffs(2); disp(num2str(c50))
            c = coeffs(3); disp(num2str(c))
            save(['D:\meso_demo\CRF_output\' animal '_curve_data_loco.mat'], "c", "c50", "rmax", "gof", "fitted_curve")

            % plot the results
            hyper_curve = figure;hold on;
            plot(x,yAvg,'r','DisplayName','test');
            plot(x,yBS,'g','DisplayName','yBS');
            legend;
            plot(fitted_curve,'b');
            str=strcat('Rmax=',num2str(coeffs(1)),' c50=',num2str(coeffs(2)),' slope=',num2str(coeffs(3)));
            annotation('textbox','String',str);
            title(animal, 'Interpreter','none')
            exportgraphics(hyper_curve,['D:\meso_demo\CRF_output\' animal '_crf_loco.pdf'],'ContentType','vector')


    
end






