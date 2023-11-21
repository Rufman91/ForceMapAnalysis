% code that travels through the folders containing colloidal-probe in the 
% name, enters, loads the experiment, goes into the folders of each 
% processed centrosome and loads the .mat file where the post-processing
% results are, plotting
% Julia Garcia Baucells 2022

format long g;
format compact;
workspace;  % make sure the workspace panel is showing
% dbstop in collect_plot_results.m at 130
% dbstop in collect_plot_results.m at 175

close all
clear
clc

path1 = '/Users/julia/MyFolder/CentrosomesExchange/Mjtraining/';
addpath('internal functions');

DirOutput = dir(path1);
FileNames = {DirOutput.name}';
NumberFolders = numel(FileNames);

msg1 = "Do you want to apply a UpperForceCutOff?";
opts1 = ["Yes" "No"];
choice1 = menu(msg1,opts1);
if choice1 == 1
    msg2 = "Which UpperForceCutOff do you want to apply?"; 
    opts2 = ["01" "02" "03"];
    choice2 = menu(msg2,opts2);
    s2 = ' ('+opts2(choice2)+')'; 
else 
    s2 = ''; 
end 

% filter out unwanted folders
pat = 'colloidal';
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
    path2 = char(strcat(path1,'/',cellstr(FileNames(ii,1)),'/','ZoomSweepFast'));
    cd(path2); % go inside date folder

    load('ZoomSweepFast.mat')
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
            load(strcat('Processed',s2,'.mat'))
            counter(n) = counter(n) + 1;
            CsEModHertz_data{n, counter(n)} = CsEModHertz(:);
            CsEModHertz_mean{n, counter(n)} = mean(CsEModHertz(:), 'omitnan');
            CsHeight_mean{n,counter(n)} = mean(CsFlatHeight(:),'omitnan');
            CsInden_mean{n, counter(n)} = mean(CsFlatInden(:),'omitnan');
            CsRadius_data{n, counter(n)} = CsRadius;
            CsAspectRatio{n,counter(n)} = mean(CsFlatHeight(:),'omitnan')/(CsRadius*2); 
            AngleThr{n,counter(n)} = T2;
            CsFlatArea_data{n, counter(n)} = CsFlatArea;
            CsVolume_data{n, counter(n)} = CsVolume; 
        end
    end
end

%% plotting data 

% % grouping of data by different days
colors = [254/255,224/255,210/255; ...
    252/255,146/255,114/255; ...
    222/255,45/255,38/255; ...
    222/255,235/255,247/255; 
    158/255,202/255,225/255; 49/255,130/255,189/255]; 
total_cs = sum(counter); 
color_map = zeros(total_cs, 3); 
for i = 1:counter(1)
    color_map(i,:) = colors(1,:); 
end 
for i = counter(1)+1:sum(counter(1:2))
    color_map(i,:) = colors(2,:); 
end 
for i = sum(counter(1:2))+1:sum(counter(1:3))
    color_map(i,:) = colors(3,:); 
end 
for i = sum(counter(1:3))+1:sum(counter(1:4))
    color_map(i,:) = colors(4,:); 
end 
for i = sum(counter(1:4))+1:sum(counter(1:5))
    color_map(i,:) = colors(5,:); 
end 
for i = sum(counter(1:5))+1:sum(counter(1:6))
    color_map(i,:) = colors(6,:); 
end 
% figure('name', 'Measurements on different days'); hold on
% colors = [0.376470588235294         0.611764705882353         0.533333333333333; ...
%     0.83921568627451         0.603921568627451         0.603921568627451; ...
%     0.650980392156863         0.650980392156863         0.650980392156863]; 
% for i = 1:length(sample_idx)
%     CsEModHertz_linear{i,1} = [];
%     BkgEModHertz_linear{i,1} = [];
%     for j = 1:counter(i)
%         plotSpread(CsEModHertz_data{i,j}.*1e-3, 'spreadWidth', 0.5, 'xNames', {'Centrosome'}, ...
%             'distributionColors', colors(i,:));
%         hold on
%         CsEModHertz_linear{i} = cat(1, CsEModHertz_linear{i,1}, CsEModHertz_data{i,j}(:));
%     end
% end
% for i = 1:length(sample_idx)
%     data_mean = mean(CsEModHertz_linear{i,1}, 'omitnan').*1e-3;
%     beeswarm(1, data_mean, 'dot_size', 5, 'MarkerEdgeColor','k', 'MarkerFaceColor', colors(i,:));
% end
% ylabel('Indentation modulus [kPa]'); box on; set(gca,'FontSize', 16, 'Linewidth', 1.5);
% ylim([0, 500]); hold off

% grouping of data by different stiffness
% total_cs = sum(counter); 
% CsEModHertz_mean_l = nan(1, total_cs); 
% for i = 1:length(sample_idx)
%     for j = 1:counter(i)
%         clr_counter = j + sum(counter(1:i-1));
%         CsEModHertz_mean_l(clr_counter) = CsEModHertz_mean{i,j}.*1e-3;
%     end 
% end 
% [B,I] = sort(CsEModHertz_mean_l); 
% color_map = zeros(total_cs, 3); 
% stiff = 4; 
% medium = 13; 
% for z = 0:(stiff-1)
%     color_map(I(end-z),:) = [0.8500 0.3250 0.0980]; % stiff
% end 
% for z = stiff:stiff+medium
%     color_map(I(end-z),:) = [0.9290 0.6940 0.1250]; % medium 
% end 
% for z = 1:sum(counter)-(stiff+medium)
%     color_map(I(z),:) = [0 0.4470 0.7410]; % soft 
% end 

figure('name', 'Volume dependence'); hold on
box on; set(gca,'FontSize', 16, 'Linewidth', 1.5);
CsEModHertz_mean_l = nan(1,total_cs);
CsVolume_l = nan(1,total_cs);
for i = 1:length(sample_idx)
    for j = 1:counter(i)
        clr_counter = j + sum(counter(1:i-1));
        c = color_map(clr_counter,:);
        CsEModHertz_mean_l(clr_counter) = CsEModHertz_mean{i,j}.*1e-3;
        CsVolume_l(clr_counter) = CsVolume_data{i,j}*1e+18;

        scatter(CsVolume_l(clr_counter), CsEModHertz_mean_l(clr_counter), 70, c, "filled");
    end
end
hold on
[xData, yData] = prepareCurveData( CsVolume_l, CsEModHertz_mean_l );

% Set up fittype and options.
ft = fittype( 'power2' );
opts = fitoptions( 'Method', 'NonlinearLeastSquares' );
opts.Display = 'Off';
opts.StartPoint = [26.8333367243072 -0.825970107131429 -62.2836917757873];

% Fit model to data.
[fitresult, gof] = fit( xData, yData, ft, opts );

% Plot fit with data.
h = plot( fitresult, xData, yData ); legend boxoff
ylabel('Indentation modulus [kPa]'); 
xlabel('Volume [\mum^3]'); xlim([-1, 10])

figure('name', 'Aspect ratio dependence'); hold on
box on; set(gca,'FontSize', 16, 'Linewidth', 1.5);
for i = 1:length(sample_idx)
    for j = 1:counter(i)
        clr_counter = j + sum(counter(1:i-1));
        c = color_map(clr_counter,:);
        scatter(CsAspectRatio{i,j}, CsEModHertz_mean{i,j}.*1e-3,  50, 'b', "filled");
        ylabel('Indentation modulus [kPa]'); 
        xlabel('Aspect ratio'); xlim([0 0.5])
    end
end

figure('name', 'Height dependence'); hold on
box on; set(gca,'FontSize', 16, 'Linewidth', 1.5);
CsEModHertz_mean_l = nan(1,total_cs);
CsHeight_mean_l = nan(1,total_cs);
for i = 1:length(sample_idx)
    for j = 1:counter(i)
        clr_counter = j + sum(counter(1:i-1));
        c = color_map(clr_counter,:);
        CsEModHertz_mean_l(clr_counter) = CsEModHertz_mean{i,j}.*1e-3;
        CsHeight_mean_l(clr_counter) = CsHeight_mean{i,j}.*1e9;

        scatter(CsHeight_mean_l(clr_counter), CsEModHertz_mean_l(clr_counter), 70, c, "filled");
    end
end

hold on
[xData, yData] = prepareCurveData( CsHeight_mean_l, CsEModHertz_mean_l );

% Set up fittype and options.
ft = fittype( 'power2' );
opts = fitoptions( 'Method', 'NonlinearLeastSquares' );
opts.Display = 'Off';
opts.StartPoint = [6133207.15536504 -2.05322611301745 -35.6361239539309];

% Fit model to data.
[fitresult, gof] = fit( xData, yData, ft, opts );

% Plot fit with data.
h = plot( fitresult, xData, yData ); legend boxoff
ylabel('Indentation modulus [kPa]'); 
xlabel('Centrosome height [nm]'); 

figure('name', 'Size dependence'); hold on
box on; set(gca,'FontSize', 16, 'Linewidth', 1.5);
for i = 1:length(sample_idx)
    for j = 1:counter(i)
        clr_counter = j + sum(counter(1:i-1));
        c = color_map(clr_counter,:);
        scatter( CsRadius_data{i,j}.*1e9, CsEModHertz_mean{i,j}.*1e-3, 70, c, "filled");
        ylabel('Indentation modulus [kPa]'); 
        xlabel('Centrosome radius [nm]'); xlim([0 2100])
    end
end

figure('name', 'Compression vs. size'); hold on
box on; set(gca,'FontSize', 16, 'Linewidth', 1.5);
for i = 1:length(sample_idx)
    for j = 1:counter(i)
        clr_counter = j + sum(counter(1:i-1));
        c = color_map(clr_counter,:);
        scatter((CsInden_mean{i,j}./CsHeight_mean{i,j})*100, CsRadius_data{i,j}.*1e9, 50, c, "filled");
        xlabel('Compression [%]')
        ylabel('Centrosome radius [nm]');
    end
end

% figure('name', 'Angle threshold'); hold on
% box on; set(gca,'FontSize', 16, 'Linewidth', 1.5);
% for i = 1:length(sample_idx)
%     for j = 1:counter(i)
%          clr_counter = j + sum(counter(1:i-1)); 
%          yyaxis left
%          data_mean_left = {AngleThr{i,j},[]};
%          plotSpread(data_mean_left,  'spreadWidth', 0.5,'xNames', {'Angle threshold', 'Flat area'}, ...
%                 'distributionMarkers',{'o','o'},'distributionColors', color_map(clr_counter,:));
%          ylabel('Angle [rad]')
%          yyaxis right
%          data_mean_right = {[],CsFlatArea_data{i,j}}; 
%                 plotSpread(data_mean_right,  'spreadWidth', 0.5,'xNames', {'Angle threshold', 'Flat area'}, ...
%                 'distributionMarkers',{'o','o'},'distributionColors', color_map(clr_counter,:));
%          ylabel('Area [m^2]')
%     end
% end

figure('name', 'Measured height and indentation depth'); hold on
box on; set(gca,'FontSize', 16, 'Linewidth', 1.5);
color_map = lines(sum(counter));
for i = 1:length(sample_idx)
    for j = 1:counter(i)
         clr_counter = j + sum(counter(1:i-1)); 
         yyaxis left
         data_mean_left = {CsInden_mean{i,j}*1e9,CsHeight_mean{i,j}*1e9, []};
         plotSpread(data_mean_left,  'spreadWidth', 0.5,'xNames', {'Indentation depth', 'Centrosome height', 'Compression'}, ...
                'distributionMarkers',{'o','o','o'},'distributionColors', color_map(clr_counter,:));
         plot([1,2], [CsInden_mean{i,j}*1e9,CsHeight_mean{i,j}*1e9], 'LineWidth',0.5, 'Color', color_map(clr_counter,:));
         ylabel('Height [nm]');

         yyaxis right
         data_mean_right = {[],[],(CsInden_mean{i,j}./CsHeight_mean{i,j})*100};
         plotSpread(data_mean_right,  'spreadWidth', 0.5,'xNames',{'Indentation depth', 'Centrosome height', 'Compression'}, ...
                'distributionMarkers',{'o','o','o'},'distributionColors', color_map(clr_counter,:));
         ylabel('% Compression'); ylim([0, 60])
    end
end
