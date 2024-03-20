% clear
% E = Experiment.load;
cd(E.ExperimentFolder)

msg1 = "Do you need to exclude force maps?";
opts1 = ["Yes" "No"];
choice1 = menu(msg1,opts1);

msg2 = "Do you want to apply a subsequent ForceMapAnalysisOptions?";
opts2 = ["Yes" "No"];
choice2 = menu(msg2,opts2);

if choice2 == 1
    msg3 = "Which ForceMapAnalysisOptions do you want to apply?";
    opts3 = ["01" "02" "03" "04"];
    choice3 = menu(msg3,opts3);
    s2 = ' ('+opts3(choice3)+')';
else
    s2 = '';
    choice3 = []; 
end

% read force map number to be excluded
if choice1 == 1
    if isfile('Exclude.txt')
        fileID = fopen('Exclude.txt', 'r');
        datacell = textscan(fileID, '%f', 'Delimiter',' ', 'CollectOutput', 1);
        fclose(fileID);
        datavalues = unique(datacell{1}); % Remove duplicates
    else
        datavalues = [];
    end
else
    datavalues = [];
end

for i = 1:E.NumForceMaps
    if ismember(i, datavalues)
        % skip evaluation round
    else
        cd(E.ForceMapFolders{i,1})
        load(strcat('Processed',s2,'.mat'))
        CsEModHertz_data{i} = CsEModHertz(:);
        CsEModHertz_mean(i) = mean(CsEModHertz(:), 'omitnan').*1e-3;
        CsFlatHeight_mean(i) = mean(CsFlatHeight(:),'omitnan').*1e9;
        CsFlatMax_data(i) = CsFlatMax*1e9; 
        CsFlatPrctile_data(i) = CsFlatPrctile*1e9; 
        CsInden_mean(i) = mean(CsFlatInden(:),'omitnan').*1e9;
        CsRadiusXY_data(i) = CsRadiusXY;
        CsAspectRatio(i) = mean(CsFlatHeight(:),'omitnan')/(CsRadiusXY*2);
        CsFlatArea_data(i) = CsFlatArea;
        CsVolume_data(i) = CsVolume*1e+18;
        CsVolumeSphereCap_data(i) = CsVolumeSphereCap*1e+18; 
    end
end

CsFlatHeight_mean(CsFlatHeight_mean == 0) = NaN;
CsFlatMax_data(CsFlatMax_data == 0) = NaN;
CsVolume_data(CsVolume_data == 0) = NaN;
CsVolumeSphereCap_data(CsVolumeSphereCap_data == 0) = NaN;
CsEModHertz_mean(CsEModHertz_mean == 0) = NaN;
nnz(~isnan(CsEModHertz_mean))

figure('name', 'Mean height dependence'); hold on
box on; set(gca,'FontSize', 16, 'Linewidth', 1.5);
for i = 1:E.NumForceMaps
    scatter(CsFlatHeight_mean(i), CsEModHertz_mean(i), 50, [0.4940 0.1840 0.5560], "filled"); %  [0 0.4470 0.7410] [0.8500 0.3250 0.0980]
end
hold on
% [xData, yData] = prepareCurveData( CsHeight_mean, CsEModHertz_mean );
% 
% % Set up fittype and options.
% ft = fittype( 'power2' );
% opts = fitoptions( 'Method', 'NonlinearLeastSquares' );
% opts.Display = 'Off';
% opts.StartPoint = [16035786.5404191 -2.23676079783161 -33.8251597804349];
% 
% % Fit model to data.
% [fitresult, gof] = fit( xData, yData, ft, opts );
% 
% % Plot fit with data.
% h = plot( fitresult, xData, yData );
ylabel('Indentation modulus [kPa]');
xlabel('Centrosome height [nm]');
% legend boxoff

figure('name', 'Maximum height dependence'); hold on
box on; set(gca,'FontSize', 16, 'Linewidth', 1.5);
for i = 1:E.NumForceMaps
    scatter(CsFlatPrctile_data(i), CsEModHertz_mean(i), 50, [0.4940 0.1840 0.5560], "filled");
end
hold on
ylabel('Indentation modulus [kPa]');
xlabel('Centrosome height [nm]');

figure('name', 'Volume dependence'); hold on
box on; set(gca,'FontSize', 16, 'Linewidth', 1.5);
for i = 1:E.NumForceMaps
    scatter(((4.*CsVolume_data(i))./(3*pi)).^(1/3)*1000, CsEModHertz_mean(i), 50, [0 0.4470 0.7410], "filled");
end
ylabel('Indentation modulus [kPa]');
xlabel('Radius [nm]');
xlim([0 1800])

figure('name', 'Maximum height vs. volume'); hold on
box on; set(gca,'FontSize', 16, 'Linewidth', 1.5);
for i = 1:E.NumForceMaps
    scatter(((4.*CsVolume_data(i))./(3*pi)).^(1/3)*1000, CsFlatPrctile_data(i), 50, [0 0.4470 0.7410], "filled");
end
ylabel('Centrosome height [nm]'); 
xlabel('Radius [nm]');

figure('name', 'Compression vs. maximum height'); hold on
box on; set(gca,'FontSize', 16, 'Linewidth', 1.5);
for i = 1:E.NumForceMaps
        scatter(CsFlatPrctile_data(i), (CsInden_mean(i)/CsFlatHeight_mean(i))*100, 50, [0.4940 0.1840 0.5560], "filled");  
end
xlabel('Centrosome height [nm]');
ylabel('Compression [%]')
