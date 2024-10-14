% clear
% E = Experiment.load;
cd(E.ExperimentFolder)
close all

msg1 = "Do you need to exclude force maps?";
opts1 = ["Yes" "No"];
choice1 = menu(msg1,opts1);

msg2 = "Do you want to apply a subsequent ForceMapAnalysisOptions?";
opts2 = ["Yes" "No"];
choice2 = menu(msg2,opts2);

if choice2 == 1
    msg3 = "Which ForceMapAnalysisOptions do you want to apply?";
    opts3 = ["01" "02" "03" "04" "05"];
    choice3 = menu(msg3,opts3);
    s2 = ' ('+opts3(choice3)+')';
else
    s2 = '';
    choice3 = []; 
end

% Read number of force maps to exclude
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
        % Skip evaluation round
    else
        cd(E.ForceMapFolders{i,1})
        load(strcat('Processed',s2,'.mat'))
        CsEModHertz_data{i} = CsEModHertz(:);
        CsEModHertz_mean(i) = mean(CsEModHertz(:), 'omitnan').*1e-3;
        CsEModHertz_std(i) = std(CsEModHertz(:), 'omitnan').*1e-3;
        CsFlatHeight_mean(i) = mean(CsFlatHeight(:),'omitnan').*1e9;
        CsFlatPrctile_data(i) = CsFlatPrctile*1e9;
        %         CsFlatMax_data(i) = CsFlatMax*1e9;
        CsInden_mean(i) = mean(CsFlatInden(:),'omitnan').*1e9;
        %         CsRadiusXY_data(i) = CsRadiusXY;
        %         CsAspectRatio(i) = mean(CsFlatHeight(:),'omitnan')/(CsRadiusXY*2);
        %         CsFlatArea_data(i) = CsFlatArea;
        CsVolume_Otsu_data(i) = CsVolume*1e+18; % From Otsu's segmentation
        %         CsVolumeSphereCap_data(i) = CsVolumeSphereCap*1e+18;

        % Calculate volume using manual segmentation
        Height{i} = E.FM{i}.get_segment_data_from_channel('Contact Height Smoothed', 'MatchString', 'Seg-02'); % Total centrosome height
        positiveHeight = max(Height{i}, 0); % Treat negative heights as zero
        Volume = sum(positiveHeight) * (E.FM{i}.ScanSizeX/E.FM{i}.NumPixelsX * E.FM{i}.ScanSizeY/E.FM{i}.NumPixelsY); % Total centrosome volume from Seg-02
        Volumes(i) = Volume*1e+18;
        %         % Equivalent radius of a sphere of the same volume
        %         EquivalentRadii_mnl(i) = ((3 * Volumes(i) / (4 * pi))^(1/3))*1000;
    end
end

CsEModHertz_mean(CsEModHertz_mean == 0) = NaN;
CsEModHertz_std(CsEModHertz_std == 0) = NaN;
CsFlatHeight_mean(CsFlatHeight_mean == 0) = NaN;
CsFlatPrctile_data(CsFlatPrctile_data == 0) = NaN; 
CsInden_mean(CsInden_mean == 0) = NaN; 
CsVolume_Otsu_data(CsVolume_Otsu_data == 0) = NaN;
Volumes(Volumes == 0) = NaN; 

% Equivalent radius of a sphere of the same volume
EquivalentRadii_mnl = (4.*Volumes./(3*pi)).^(1/3)*1000; 
EquivalentRadii_auto = (4.*CsVolume_Otsu_data./(3*pi)).^(1/3)*1000; 

% figure('name', 'Maximum height dependence'); hold on
% box on; set(gca,'FontSize', 16, 'Linewidth', 1.5);
% for i = 1:E.NumForceMaps
%     scatter(CsFlatPrctile_data(i), CsEModHertz_mean(i), 50, [0.4940 0.1840 0.5560], "filled");
% end
% hold on
% ylabel('Indentation modulus [kPa]');
% xlabel('Centrosome height [nm]');

if isempty(choice3)
    c =  [55/255 126/255 184/255];
elseif choice3 == 1
    c =  [152/255 78/255 163/255]; % Thin film (not bonded)
elseif choice3 == 2
    c = [255/255 127/255 0/255]; % Topography
elseif choice3 == 3
    c = [77/255 175/255 74/255]; % Thin film (not bonded) + topography
elseif choice3 == 4
    c = [228/255 26/255 28/255]; % Thin film (bonded)
elseif choice3 == 5
    c = [153/255 153/255 153/255]; % Thin film (bonded) + topography
end

figure('name', 'Centrosome volume dependence'); hold on
box on; set(gca,'FontSize', 18, 'Linewidth', 1.5);
scatter(EquivalentRadii_mnl, CsEModHertz_mean, 60, c, "filled");
errorbar(EquivalentRadii_mnl, CsEModHertz_mean, CsEModHertz_std, 'o', 'Color', c);
ylabel('Indentation modulus [kPa]');
xlabel('Centrosome equivalent radius [nm]');
xlim([0 1600]); ylim([-50 350])

RobustFit = false;

% Find indices where any of the arrays have NaN values
nanIndices = isnan(EquivalentRadii_mnl) | isnan(CsEModHertz_mean);

% Remove NaN values from each array
EquivalentRadii_clean = EquivalentRadii_mnl(~nanIndices);
CsEModHertz_mean_clean = CsEModHertz_mean(~nanIndices);
addCorrelationInfo(EquivalentRadii_clean, CsEModHertz_mean_clean, RobustFit);

figure('name', 'Centrosome height dependence'); hold on
box on; set(gca,'FontSize', 18, 'Linewidth', 1.5);
scatter(CsFlatPrctile_data, CsEModHertz_mean, 60, c, "filled");
errorbar(CsFlatPrctile_data, CsEModHertz_mean, CsEModHertz_std, 'o', 'Color', c);

% hold on
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

% % Plot fit with data.
% h = plot( fitresult, xData, yData );
ylabel('Indentation modulus [kPa]');
xlabel('Centrosome max. height [nm]');
ylim([-50 350]); xlim([0, 1100])

RobustFit = false;

% Find indices where any of the arrays have NaN values
nanIndices = isnan(CsFlatPrctile_data) | isnan(CsEModHertz_mean);

% Remove NaN values from each array
CsFlatPrctile_data_clean = CsFlatPrctile_data(~nanIndices);
CsEModHertz_mean_clean = CsEModHertz_mean(~nanIndices);
addCorrelationInfo(CsFlatPrctile_data_clean, CsEModHertz_mean_clean, RobustFit);

figure('name', 'Centrosome maximum height vs. Manual volume'); hold on
box on; set(gca,'FontSize', 18, 'Linewidth', 1.5);
scatter(EquivalentRadii_mnl, CsFlatPrctile_data, 60, [0 0.4470 0.7410], "filled");
ylabel('Centrosome maximum height [nm]'); 
xlabel('Centrosome equivalent radius [nm] - Manual segmentation');
ylim([0, 1100]); xlim([0 1800])

% Linear regression fit
validIdx = ~isnan(EquivalentRadii_mnl) & ~isnan(CsFlatPrctile_data); % Remove NaN values from both variables
EquivalentRadii_valid = EquivalentRadii_mnl(validIdx);
CsFlatPrctile_data_valid = CsFlatPrctile_data(validIdx);
[p, S] = polyfit(EquivalentRadii_valid, CsFlatPrctile_data_valid, 1);  % Linear fit
yfit = polyval(p, EquivalentRadii_valid);  % Predicted values based on fit
correlation_coefficient = corr(EquivalentRadii_valid', CsFlatPrctile_data_valid', 'Type', 'Pearson');

% Calculate the confidence intervals for the slope and intercept
[fit_ci, ~] = polyconf(p, EquivalentRadii_valid, S, 'alpha', 0.05); % 95% CI for the fit

% Calculate R-squared
SS_res = sum((CsFlatPrctile_data_valid - yfit).^2);  % Sum of squares of residuals
SS_tot = sum((CsFlatPrctile_data_valid - mean(CsFlatPrctile_data_valid)).^2);  % Total sum of squares
R_squared = 1 - (SS_res / SS_tot);

% Linear equation
slope = p(1);
intercept = p(2);

% Add the linear regression line
hold on;
plot(EquivalentRadii_valid, yfit, '-r', 'LineWidth', 1);  % Plot fit line
% Plot the confidence intervals for the fit line
plot(EquivalentRadii_valid, fit_ci, 'r--', 'LineWidth', 1);  % Plot confidence intervals
% Display correlation, R-squared, and linear equation in the legend
legend('Data points', ['Linear fit: y = ' num2str(slope, '%.2f') 'x + ' num2str(intercept, '%.2f')], ['R^2 = ', num2str(R_squared, '%.2f')]); legend boxoff
hold off

figure('name', 'Centrosome maximum height vs. Otsu volume'); hold on
box on; set(gca,'FontSize', 18, 'Linewidth', 1.5);
scatter(EquivalentRadii_auto, CsFlatPrctile_data, 60, [0 0.4470 0.7410], "filled");
ylabel('Centrosome maximum height [nm]'); 
xlabel('Centrosome equivalent radius [nm] - Otsu thresholding');
ylim([0, 1100]); xlim([0 1800])

% Linear regression fit
validIdx = ~isnan(EquivalentRadii_auto) & ~isnan(CsFlatPrctile_data); % Remove NaN values from both variables
EquivalentRadii_valid = EquivalentRadii_auto(validIdx);
CsFlatPrctile_data_valid = CsFlatPrctile_data(validIdx);
[p, S] = polyfit(EquivalentRadii_valid, CsFlatPrctile_data_valid, 1);  % Linear fit
yfit = polyval(p, EquivalentRadii_valid);  % Predicted values based on fit
correlation_coefficient = corr(EquivalentRadii_valid', CsFlatPrctile_data_valid', 'Type', 'Pearson');

% Calculate the confidence intervals for the slope and intercept
[fit_ci, delta] = polyconf(p, EquivalentRadii_valid, S, 'alpha', 0.05); % 95% CI for the fit

% Calculate R-squared
SS_res = sum((CsFlatPrctile_data_valid - yfit).^2);  % Sum of squares of residuals
SS_tot = sum((CsFlatPrctile_data_valid - mean(CsFlatPrctile_data_valid)).^2);  % Total sum of squares
R_squared = 1 - (SS_res / SS_tot);

% Linear equation
slope = p(1);
intercept = p(2);

% Add the linear regression line
hold on;
plot(EquivalentRadii_valid, yfit, '-r', 'LineWidth', 1);  % Plot fit line
% Plot the confidence intervals for the fit line
plot(EquivalentRadii_valid, fit_ci, 'r--', 'LineWidth', 1);  % Plot confidence intervals
% Display correlation, R-squared, and linear equation in the legend
legend('Data points', ['Linear fit: y = ' num2str(slope, '%.2f') 'x + ' num2str(intercept, '%.2f')], ['R^2 = ', num2str(R_squared, '%.2f')]); legend boxoff
hold off

 
figure('name', 'Compression vs. height'); hold on
box on; set(gca,'FontSize', 16, 'Linewidth', 1.5);
for i = 1:E.NumForceMaps
        Compression = (CsInden_mean(i)/CsFlatHeight_mean(i))*100; 
        scatter(CsFlatPrctile_data(i), Compression, 60, c, "filled");  
end
xlabel('Centrosome maximum height [nm]');
ylabel('Compression [%]')
xlim([0, 1100]); ylim([0 70])


figure('name', 'Compression vs. height'); 
hold on;
box on; 
set(gca,'FontSize', 16, 'Linewidth', 1.5);

lowCompIdx = [];
midCompIdx = [];
highCompIdx = [];
edgeColor = [153/255 153/255 153/255];
lowColor = [145/255 207/255 96/255];
midColor = [255/255 255/255 191/255];
highColor = [252/255 141/255 89/255];

for i = 1:E.NumForceMaps
    Compression = (CsInden_mean(i) / CsFlatHeight_mean(i)) * 100;
    
    % Color code based on Compression value
    if Compression < 25
        c = lowColor; % Green for Compression < 25%
        lowCompIdx = [lowCompIdx i]; 
    elseif Compression >= 25 && Compression <= 35
        c = midColor; % Yellow for Compression between 25-35%
        midCompIdx = [midCompIdx i]; 
    else
        c = highColor; % Red for Compression > 35%
        highCompIdx = [highCompIdx i];
    end
    
    % Plot the data point
  scatter(CsFlatPrctile_data(i), CsEModHertz_mean(i), 60, 'MarkerEdgeColor', edgeColor, 'MarkerFaceColor', c);
  errorbar(CsFlatPrctile_data, CsEModHertz_mean, CsEModHertz_std, 'o', 'Color', edgeColor);
end

ylabel('Indentation modulus [kPa]');
xlabel('Centrosome max. height [nm]');
ylim([-50 350]); xlim([0, 1100])
% Add dummy scatter plots for legend
h1 = scatter(nan, nan, 60, 'MarkerEdgeColor', edgeColor, 'MarkerFaceColor', lowColor);
h2 = scatter(nan, nan, 60, 'MarkerEdgeColor', edgeColor, 'MarkerFaceColor', midColor);
h3 = scatter(nan, nan, 60, 'MarkerEdgeColor', edgeColor, 'MarkerFaceColor', highColor);
% Create legend
legend([h1, h2, h3], '< 25% Compression', '25-35% Compression', '> 35% Compression', 'Location', 'best');


function addCorrelationInfo(x, y, useRobustFit)
    % Calculate correlation coefficient and p-value
    [r, p] = corrcoef(x, y);
    r_value = r(1, 2);
    p_value = p(1, 2);
    
    % Calculate the number of samples needed for 80% power at alpha = 0.05
    Z_alpha_over_2 = 1.96; % for alpha = 0.05
    Z_beta = 0.84; % for 80% power
    fisher_z = 0.5 * log((1 + r_value) / (1 - r_value));
    needed_samples = (Z_alpha_over_2 + Z_beta)^2 / fisher_z^2 + 3;
    
    % Display correlation line
    hold on;
    if useRobustFit
        coeffs = robustfit(x, y);
        fittedX = linspace(min(x), max(x), 200);
        fittedY = coeffs(1) + coeffs(2) * fittedX;
    else
        coeffs = polyfit(x, y, 1);
        fittedX = linspace(min(x), max(x), 200);
        fittedY = polyval(coeffs, fittedX);
    end
    plot(fittedX, fittedY, 'r-', 'LineWidth', 1);
    
    % Display r, p values, and needed samples
    text(min(x) + 0.05 * (max(x) - min(x)), max(y) - 0.05 * (max(y) - min(y)), ...
         sprintf('r = %.2f\np = %.2g\nNeeded samples = %.0f', r_value, p_value, ceil(needed_samples)), ...
         'EdgeColor', 'black', 'BackgroundColor', 'white', 'Margin', 5);
    hold off;
end
