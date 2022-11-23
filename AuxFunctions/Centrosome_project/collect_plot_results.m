% code that travels through the folders containing colloidal-probe in the 
% name, enters, loads the experiment, goes into the folders of each 
% processed centrosome and loads the .mat file where the post-processing
% results are, plots the results in a spread plot
% Julia Garcia Baucells 2022

format long g;
format compact;
workspace;  % make sure the workspace panel is showing
% dbstop in collect_plot_results.m at 123
dbstop in collect_plot_results.m at 162
% dbstop in collect_plot_results.m at 184
% dbstop in collect_plot_results.m at 213

close all
clear
clc

path1 = '/Users/julia/MyFolder/CentrosomesExchange/Mjtraining/';
addpath('internal functions');

DirOutput = dir(path1);
FileNames = {DirOutput.name}';
NumberFolders = numel(FileNames);

% filter out unwanted folders
pat = 'colloidal-probe';
counter = 1;
sample_idx = nan(NumberFolders,1);
for m = 1:NumberFolders
    TF = contains(FileNames(m,1),pat);
    if TF == 1
        sample_idx(counter) = m;
    else
        sample_idx(counter) = 0;
    end
    counter = counter + 1;
end
sample_idx = nonzeros(sample_idx);
clear TF; clear counter;

counter(1:length(sample_idx)) = 0; 
for n = 1:length(sample_idx)
    ii = sample_idx(n);
    path2 = char(strcat(path1,'/',cellstr(FileNames(ii,1)),'/','E'));
    cd(path2); % go inside date folder

    load('E.mat')
    if isfile('Exclude.txt')
        fileID = fopen('Exclude.txt', 'r');
        datacell = textscan(fileID, '%f', 'Delimiter',' ', 'CollectOutput', 1);
        fclose(fileID);
        datavalues = datacell{1};    %as a numeric array
    end
   
    % collect data
    for i = 1:obj.NumForceMaps
        if ismember(i, datavalues)
            % skip evaluation round
        else
            cd(obj.ForceMapFolders{i,1})
            load('Processed.mat')
            counter(n) = counter(n) + 1;
            CsEModHertz_data{n, counter(n)} = CsEModHertz(:);
            BkgEModHertz_data{n, counter(n)} = BkgEModHertz(:);
            CsEModHertz_mean{n, counter(n)} = nanmean(CsEModHertz(:));
            CsHeight_mean{n,counter(n)} = nanmean(CsFlatHeight(:));
            CsInden_mean{n, counter(n)} = nanmean(CsFlatInden(:));
            CsRadius_data{n, counter(n)} = CsRadius;
            CsAspectRatio{n,counter(n)} = nanmean((CsRadius*2)./CsFlatHeight(:)); % nanmean(CsFlatHeight(:))./(CsRadius*2)
            AngleThr{n,counter(n)} = T2;
            CsFlatArea_data{n, counter(n)} = CsFlatArea;
        end
    end
end

% plot data
figure('name', 'Measurements on different days'); hold on
colors = [0.376470588235294         0.611764705882353         0.533333333333333; ...
    0.83921568627451         0.603921568627451         0.603921568627451; ...
    0.650980392156863         0.650980392156863         0.650980392156863]; 
for i = 1:length(sample_idx)
    CsEModHertz_linear{i,1} = [];
    BkgEModHertz_linear{i,1} = [];
    for j = 1:counter(i)
        plotSpread({CsEModHertz_data{i,j},BkgEModHertz_data{i,j}} , 'spreadWidth', 0.5, 'xNames', {'Centrosome', 'Background'}, ...
            'distributionColors', colors(i,:));
        hold on
        CsEModHertz_linear{i} = cat(1, CsEModHertz_linear{i,1}, CsEModHertz_data{i,j}(:));
        BkgEModHertz_linear{i} = cat(1, BkgEModHertz_linear{i,1}, BkgEModHertz_data{i,j}(:));
    end
end
for i = 1:length(sample_idx)
    data_mean = cat(2,nanmean(CsEModHertz_linear{i,1}), nanmean(BkgEModHertz_linear{i,1}));
    beeswarm([1,2], data_mean, 'dot_size', 5, 'MarkerEdgeColor','k', 'MarkerFaceColor', colors(i,:));
end
ylabel('Indentation modulus [Pa]'); box on; set(gca,'FontSize', 16, 'Linewidth', 1.5);
ylim([0, 4.5*1e6])

figure('name', 'Mechanics geometrical dependence'); hold on
box on; set(gca,'FontSize', 16, 'Linewidth', 1.5);
color_map = hot(sum(counter));
for i = 1:length(sample_idx)
    for j = 1:counter(i)
        plotSpread({CsEModHertz_data{i,j},[]} , 'spreadWidth', 0.5, 'xNames', {'Centrosome stiffness', 'Centrosome aspect ratio'}, ...
            'distributionColors', [0.5 0.5 0.5]);
        hold on
    end
    for j = 1:counter(i)
        clr_counter = j + sum(counter(1:i-1));
        CsAspectRatio_l(clr_counter) = CsAspectRatio{i,j};
        CsEModHertz_mean_l(clr_counter) = CsEModHertz_mean{i,j};

        yyaxis left
        data_mean_left = cat(2, CsEModHertz_mean{i,j},NaN);
        beeswarm([1,2], data_mean_left, 'dot_size', 5, 'MarkerEdgeColor','k', 'MarkerFaceColor', color_map(clr_counter,:));
        ylabel('Indentation modulus [Pa]')
        ylim([0, 3*1e5])

        yyaxis right
        data_mean_right = cat(2, NaN, CsAspectRatio{i,j});
        beeswarm([1,2], data_mean_right, 'dot_size', 5, 'MarkerEdgeColor','k', 'MarkerFaceColor', color_map(clr_counter,:));
        ylabel('Aspect ratio');
        ylim([0,8]) % ylim([0,0.6]) 
    end
end
% paired-sample t-test
[h1,p2] = ttest2(CsAspectRatio_l,CsEModHertz_mean_l); % mean data


[B,I] = sort(CsEModHertz_mean_l); 
% color_mapv1 = summer(sum(counter));
color_map = zeros(sum(counter), 3); 
stiff = 5; 
medium = 3; 
for z = 0:(stiff-1)
    color_map(I(end-z),:) = [0.8500 0.3250 0.0980]; % stiff
end 
clear z 
for z = stiff:stiff+medium
    color_map(I(end-z),:) = [0.9290 0.6940 0.1250]; % medium 
end 
for z = 1:sum(counter)-(stiff+medium)
    color_map(I(z),:) = [0 0.4470 0.7410]; % soft 
end 

figure('name', 'Mechanics height dependence'); hold on
box on; set(gca,'FontSize', 16, 'Linewidth', 1.5);
for i = 1:length(sample_idx)
    for j = 1:counter(i)
        clr_counter = j + sum(counter(1:i-1));
        yyaxis left
        plotSpread({CsEModHertz_data{i,j}.*1e-3,[]} , 'spreadWidth', 0.5, 'xNames', {'Centrosome stiffness', 'Centrosome height'}, ...
            'distributionColors', color_map(clr_counter,:));
        hold on

        CsHeight_mean_l(clr_counter) = CsHeight_mean{i,j};
        CsEModHertz_mean_l(clr_counter) = CsEModHertz_mean{i,j};


        data_mean_left = cat(2, CsEModHertz_mean{i,j}.*1e-3,NaN);
        beeswarm([1,2], data_mean_left, 'dot_size', 5, 'MarkerEdgeColor','k', 'MarkerFaceColor', color_map(clr_counter,:));
        ylabel('Indentation modulus [kPa]')
        ylim([0, 3*1e2]) % 3*1e5

        yyaxis right
        data_mean_right = cat(2, NaN, CsHeight_mean{i,j}.*1e9);
        beeswarm([1,2], data_mean_right, 'dot_size', 5, 'MarkerEdgeColor','k', 'MarkerFaceColor', color_map(clr_counter,:));
        ylabel('Centrosome height [nm]');
        ylim([0,800])
    end
end
% paired-sample t-test
[h2,p2] = ttest2(CsHeight_mean_l,CsEModHertz_mean_l); % mean data

figure('name', 'Mechanics size dependence'); hold on
box on; set(gca,'FontSize', 16, 'Linewidth', 1.5);
color_map = hot(sum(counter));
for i = 1:length(sample_idx)
    for j = 1:counter(i)
        plotSpread({CsEModHertz_data{i,j},[]} , 'spreadWidth', 0.5, 'xNames', {'Centrosome stiffness', 'Centrosome size'}, ...
            'distributionColors', [0.5 0.5 0.5]);
        hold on
    end
    for j = 1:counter(i)
        clr_counter = j + sum(counter(1:i-1));
        CsRadius_data_l(clr_counter) = CsRadius_data{i,j};
        CsEModHertz_mean_l(clr_counter) = CsEModHertz_mean{i,j};

        yyaxis left
        data_mean_left = cat(2, CsEModHertz_mean{i,j},NaN);
        beeswarm([1,2], data_mean_left, 'dot_size', 5, 'MarkerEdgeColor','k', 'MarkerFaceColor', color_map(clr_counter,:));
        ylabel('Indentation modulus [Pa]')
        ylim([0, 3*1e5])

        yyaxis right
        data_mean_right = cat(2, NaN, CsRadius_data{i,j});
        beeswarm([1,2], data_mean_right, 'dot_size', 5, 'MarkerEdgeColor','k', 'MarkerFaceColor', color_map(clr_counter,:));
        ylabel('Centrosome radius [m]');
        ylim([0,1.5e-6])
    end
end
% paired-sample t-test
[h3,p3] = ttest2(CsRadius_data_l,CsEModHertz_mean_l); % mean data

figure('name', 'Compression size dependence'); hold on
box on; set(gca,'FontSize', 16, 'Linewidth', 1.5);
color_map = hot(sum(counter));
for i = 1:length(sample_idx)
     plotSpread({[],[]} , 'spreadWidth', 0.5, 'xNames', {'Centrosome size', 'Compression'}, ...
            'distributionColors', [0.5 0.5 0.5]);
        hold on
    for j = 1:counter(i)
        clr_counter = j + sum(counter(1:i-1));
        Compression_data_l(clr_counter) = (CsInden_mean{i,j}./CsHeight_mean{i,j})*100;
        CsRadius_data_l(clr_counter) = CsRadius_data{i,j};

        yyaxis left
        data_mean_left = cat(2, CsRadius_data{i,j},NaN);
        beeswarm([1,2], data_mean_left, 'dot_size', 5, 'MarkerEdgeColor','k', 'MarkerFaceColor', color_map(clr_counter,:));
        ylabel('Centrosome radius [m]');
        ylim([0,1.5e-6])

        yyaxis right
        data_mean_right = cat(2, NaN, (CsInden_mean{i,j}./CsHeight_mean{i,j})*100);
        beeswarm([1,2], data_mean_right, 'dot_size', 5, 'MarkerEdgeColor','k', 'MarkerFaceColor', color_map(clr_counter,:));
        ylabel('Compression [%]');
        ylim([20,60])
    end
end
% paired-sample t-test
[h4,p4] = ttest2(Compression_data_l,CsRadius_data_l); % mean data

figure('name', 'Angle threshold'); hold on
box on; set(gca,'FontSize', 16, 'Linewidth', 1.5);
for i = 1:length(sample_idx)
    for j = 1:counter(i)
         clr_counter = j + sum(counter(1:i-1)); 
         yyaxis left
         data_mean_left = {AngleThr{i,j},[]};
         plotSpread(data_mean_left,  'spreadWidth', 0.5,'xNames', {'Angle threshold', 'Flat area'}, ...
                'distributionMarkers',{'o','o'},'distributionColors', color_map(clr_counter,:));
         ylabel('Angle [rad]')
         yyaxis right
         data_mean_right = {[],CsFlatArea_data{i,j}}; 
                plotSpread(data_mean_right,  'spreadWidth', 0.5,'xNames', {'Angle threshold', 'Flat area'}, ...
                'distributionMarkers',{'o','o'},'distributionColors', color_map(clr_counter,:));
         ylabel('Area [m^2]')
    end
end

figure('name', 'Measured height and indentation depth'); hold on
box on; set(gca,'FontSize', 16, 'Linewidth', 1.5);
color_map = summer(sum(counter));
for i = 1:length(sample_idx)
    for j = 1:counter(i)
         clr_counter = j + sum(counter(1:i-1)); 
         yyaxis left
         data_mean_left = {CsInden_mean{i,j},CsHeight_mean{i,j}, []};
         plotSpread(data_mean_left,  'spreadWidth', 0.5,'xNames', {'Indentation depth', 'Centrosome height', 'Compression'}, ...
                'distributionMarkers',{'o','o','o'},'distributionColors', color_map(clr_counter,:));
         plot([1,2], [CsInden_mean{i,j},CsHeight_mean{i,j}], 'LineWidth',0.5, 'Color', color_map(clr_counter,:));
         ylabel('Height [m]');

         yyaxis right
         data_mean_right = {[],[],(CsInden_mean{i,j}./CsHeight_mean{i,j})*100};
         plotSpread(data_mean_right,  'spreadWidth', 0.5,'xNames',{'Indentation depth', 'Centrosome height', 'Compression'}, ...
                'distributionMarkers',{'o','o','o'},'distributionColors', color_map(clr_counter,:));
         ylabel('% Compression');
    end
end
