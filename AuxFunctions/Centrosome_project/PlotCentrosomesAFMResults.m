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
        CsInden_std(i) = std(CsFlatInden(:),'omitnan').*1e9;
%         CsEffectiveRadius_mean(i) = mean(CsEffectiveRadius(:),'omitnan').*1e9;
%         CsEffectiveRadius_std(i) = std(CsEffectiveRadius(:),'omitnan').*1e9;
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
CsInden_std(CsInden_std == 0) = NaN; 
% CsEffectiveRadius_mean(CsEffectiveRadius_mean == 0) = NaN; 
% CsEffectiveRadius_std(CsEffectiveRadius_std == 0) = NaN; 
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
    c =  [255/255 127/255 0/255];
elseif choice3 == 1
    c =  [152/255 78/255 163/255]; % Thin film (not bonded)
elseif choice3 == 2
    c = [55/255 126/255 184/255]; % Topography
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
xlim([0 1500]); ylim([0 400])

%%%% Color-code based on centrosome radius 
figure('Name', 'Centrosome volume dependence'); 
hold on;
box on; 
set(gca,'FontSize', 18, 'Linewidth', 1.5);

% Create color vector based on radius threshold
clrs = zeros(length(EquivalentRadii_mnl), 3); % Initialize color matrix
orange = [0.8500 0.3250 0.0980]; % MATLAB default orange
blue = [0 0.4470 0.7410]; % MATLAB default blue

filtered_radii = [];
filtered_moduli = [];
filtered_stds = [];
filtered_colors = [];

% Assign colors based on radius
for i = 1:length(EquivalentRadii_mnl)
    if EquivalentRadii_mnl(i) > 500
        clrs(i,:) = orange;
    else
        clrs(i,:) = blue;
    end
end

% Add error bars with matching colors
for i = 1:length(EquivalentRadii_mnl)
    % Skip indices 14-32 (Day 2)
    if i >= 14 && i <= 32
        continue;
    end
    % Add to filtered arrays
    filtered_radii(end+1) = EquivalentRadii_mnl(i);
    filtered_moduli(end+1) = CsEModHertz_mean(i);
    filtered_stds(end+1) = CsEModHertz_std(i);

    % Assign color based on radius
    if EquivalentRadii_mnl(i) > 500
        filtered_colors(end+1,:) = orange;
    else
        filtered_colors(end+1,:) = blue;
    end
end

% Create scatter plot with filtered data
scatter(filtered_radii, filtered_moduli, 60, filtered_colors, "filled");

% Add error bars with matching colors
for i = 1:length(filtered_radii)
    if filtered_radii(i) > 500
        errorbar(filtered_radii(i), filtered_moduli(i), filtered_stds(i), ...
                'o', 'Color', orange, 'MarkerFaceColor', orange);
    else
        errorbar(filtered_radii(i), filtered_moduli(i), filtered_stds(i), ...
                'o', 'Color', blue, 'MarkerFaceColor', blue);
    end
end

% Add labels and limits
ylabel('Indentation modulus [kPa]');
xlabel('Centrosome equivalent radius [nm]');
xlim([0 1500]); 
ylim([0 350]);

% Add legend
h = zeros(2,1);
h(1) = plot(NaN,NaN,'o','MarkerEdgeColor',blue,'MarkerFaceColor',blue);
h(2) = plot(NaN,NaN,'o','MarkerEdgeColor',orange,'MarkerFaceColor',orange);
legend(h, {'Radius ≤ 500 nm', 'Radius > 500 nm'}, 'Location', 'best');
legend box off

% % Optional: Add reference line at 500 nm
% xline(500, '--k', 'LineWidth', 1, 'Alpha', 0.5);

% Count non-NaN values in original data
original_nonnan = sum(~isnan(EquivalentRadii_mnl) & ~isnan(CsEModHertz_mean) & ~isnan(CsEModHertz_std));
fprintf('Original data points (non-NaN): %d\n', original_nonnan);

% Count non-NaN values in filtered data
filtered_nonnan = sum(~isnan(filtered_radii) & ~isnan(filtered_moduli) & ~isnan(filtered_stds));
fprintf('Filtered data points (non-NaN): %d\n', filtered_nonnan);

% Count how many were excluded by index filter
excluded_by_index = sum((1:length(EquivalentRadii_mnl)) >= 14 & (1:length(EquivalentRadii_mnl)) <= 32);
fprintf('Points excluded by index filter: %d\n', excluded_by_index);

% Count how many non-NaN points were excluded
valid_indices = ~isnan(EquivalentRadii_mnl) & ~isnan(CsEModHertz_mean) & ~isnan(CsEModHertz_std);
excluded_nonnan = sum(valid_indices & ((1:length(EquivalentRadii_mnl)) >= 14 & (1:length(EquivalentRadii_mnl)) <= 32));
fprintf('Valid non-NaN points excluded: %d\n', excluded_nonnan);

%%%%

RobustFit = false;

% Find indices where any of the arrays have NaN values
nanIndices = isnan(filtered_radii) | isnan(filtered_moduli);

% Remove NaN values from each array
EquivalentRadii_clean = filtered_radii(~nanIndices);
CsEModHertz_mean_clean = filtered_moduli(~nanIndices);

% Calculate distance correlation
x = EquivalentRadii_clean'; 
y = CsEModHertz_mean_clean'; 
dcor_observed = distcorr(x, y);

% Permutation test 
num_permutations = 10000; 
dcor_permuted = zeros(1, num_permutations);

% Permutation test
for i = 1:num_permutations
    % Randomly shuffle the y-values
    y_permuted = y(randperm(length(y)));
    
    % Calculate the distance correlation for the permuted data
    dcor_permuted(i) = distcorr(x, y_permuted);
end

% Calculate the p-value
p_value = mean(dcor_permuted >= dcor_observed);

% Display the results
fprintf('Observed distance correlation: %.4f\n', dcor_observed);
fprintf('p-value: %.4f\n', p_value);

% Calculate linear correlation
addCorrelationInfo(EquivalentRadii_clean, CsEModHertz_mean_clean, RobustFit);

% Get indices for each group
idx_large = filtered_radii > 500;
idx_small = filtered_radii <= 500;

% Extract moduli for each group
moduli_large = filtered_moduli(idx_large);
moduli_small = filtered_moduli(idx_small);

mean_large = mean(moduli_large, 'omitnan');
mean_small = mean(moduli_small, 'omitnan');

std_large = std(moduli_large, 'omitnan');
std_small = std(moduli_small, 'omitnan');

fold_change = mean_small / mean_large;

% Check normality
[h_large, p_large] = lillietest(moduli_large);
[h_small, p_small] = lillietest(moduli_small);

if p_large >= 0.05 && p_small >= 0.05
    % Both normal → check variances
    [~, p_var] = vartest2(moduli_large, moduli_small);
    if p_var >= 0.05
        % Equal variances → standard t-test
        [h, p] = ttest2(moduli_large, moduli_small);
        test_used = 'Student''s t-test (equal variances)';
    else
        % Unequal variances → Welch’s t-test
        [h, p] = ttest2(moduli_large, moduli_small, 'Vartype', 'unequal');
        test_used = 'Welch''s t-test (unequal variances)';
    end
else
    % Non-normal → Mann-Whitney U test
    [p, h] = ranksum(moduli_large, moduli_small);
    test_used = 'Mann-Whitney U test';
end

fprintf('\nSelected test: %s\n', test_used);
fprintf('p-value = %.4f\n', p);
if h == 1
    fprintf('Conclusion: Significant difference (p < 0.05)\n');
else
    fprintf('Conclusion: No significant difference (p ≥ 0.05)\n');
end

% figure('name', 'Centrosome height dependence'); hold on
% box on; set(gca,'FontSize', 18, 'Linewidth', 1.5);
% scatter(CsFlatPrctile_data, CsEModHertz_mean, 60, c, "filled");
% errorbar(CsFlatPrctile_data, CsEModHertz_mean, CsEModHertz_std, 'o', 'Color', c);

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
% ylabel('Indentation modulus [kPa]');
% xlabel('Centrosome max. height [nm]');
% ylim([-50 350]); xlim([0, 1100])

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

figure('name', 'Compression Indentation modulus dependence'); hold on
box on; set(gca,'FontSize', 16, 'Linewidth', 1.5);
for i = 1:E.NumForceMaps
    % Skip indices 14-32 (Day 2)
    if i >= 14 && i <= 32
        continue;
    end
    Compression = (CsInden_mean(i)/CsFlatHeight_mean(i))*100;
    scatter(EquivalentRadii_mnl(i), Compression, 60, c, "filled");
end
xlabel('Centrosome equivalent radius [nm]');
ylabel('Compression [%]')
xlim([0, 1500]); 

%% Color-code based on compression - IMPROVED VERSION
figure('name', 'Compression vs. height'); 
hold on; 
box on; 
set(gca,'FontSize', 18, 'Linewidth', 1.5); % Increased font size to match previous plot

% Define colors (using more vibrant versions of your colors)
lowColor  = [0.18, 0.55, 0.34]; 
midColor  = [0.93, 0.69, 0.13]; 
highColor = [0.80, 0.25, 0.15]; 

% Create logical index of points to include (not between 14-32)
includeIdx = true(1, E.NumForceMaps);
includeIdx(14:32) = false;

% Initialize arrays for plotting
plot_x = [];
plot_y = [];
plot_err = [];
plot_colors = [];

for i = 1:E.NumForceMaps
    % Skip indices 14-32 (Day 2)
    if ~includeIdx(i)
        continue;
    end

    Compression = (CsInden_mean(i) / CsFlatHeight_mean(i)) * 100;

    % Color code based on Compression value
    if Compression < 25
        c = lowColor;
    elseif Compression >= 25 && Compression <= 35
        c = midColor;
    else
        c = highColor;
    end
    
    % Store values for plotting
    plot_x = [plot_x EquivalentRadii_mnl(i)];
    plot_y = [plot_y CsEModHertz_mean(i)];
    plot_err = [plot_err CsEModHertz_std(i)];
    plot_colors = [plot_colors; c];
end

% Create scatter plot with filled markers (no edge)
scatter(plot_x, plot_y, 60, plot_colors, "filled");

% Add error bars with matching colors (no marker)
for i = 1:length(plot_x)
    errorbar(plot_x(i), plot_y(i), plot_err(i), 'o', ...
             'Color', plot_colors(i,:), 'MarkerFaceColor', plot_colors(i,:));
end

ylabel('Indentation modulus [kPa]');
xlabel('Centrosome equivalent radius [nm]'); % Changed to match previous plot
ylim([0 350]); xlim([0 1500]); % Consistent limits with previous plot

% Create legend with custom markers
h = zeros(3,1);
h(1) = plot(NaN,NaN,'o','MarkerEdgeColor',lowColor,'MarkerFaceColor',lowColor,'MarkerSize',8);
h(2) = plot(NaN,NaN,'o','MarkerEdgeColor',midColor,'MarkerFaceColor',midColor,'MarkerSize',8);
h(3) = plot(NaN,NaN,'o','MarkerEdgeColor',highColor,'MarkerFaceColor',highColor,'MarkerSize',8);
legend(h, {'< 25% Compression', '25-35% Compression', '> 35% Compression'}, ...
       'Location', 'best', 'FontSize', 16);
legend box off;

%% Color-code based on acquisition day
figure('name', 'Acquisition day/cantilever tip dependence'); 
hold on; box on; 
set(gca,'FontSize', 16, 'Linewidth', 1.5);

day1Color = [224/255 236/255 244/255];
day2Color = [158/255 188/255 218/255];
day3Color = [136/255 86/255 167/255];
day4Color = [254/255 232/255 200/255];
day5Color = [253/255 187/255 132/255];
day6Color = [227/255 74/255 51/255];
for i = 1:E.NumForceMaps
    % Color code based on acquisition day
    if i >= 1 && i <= 13
        c = day1Color;
    elseif i >= 14 && i <= 32
        c = day2Color;
    elseif i >= 33 && i <= 45
        c = day3Color;
    elseif i >= 46 && i <= 47
        c = day4Color;
    elseif i >= 48 && i <= 52
        c = day5Color;
    elseif i >= 53
        c = day6Color;
    end

    % Plot the data point
    scatter(CsFlatPrctile_data(i), CsEModHertz_mean(i), 60, 'MarkerEdgeColor', edgeColor, 'MarkerFaceColor', c);
    errorbar(CsFlatPrctile_data, CsEModHertz_mean, CsEModHertz_std, 'o', 'Color', edgeColor);
end

ylabel('Indentation modulus [kPa]');
xlabel('Centrosome max. height [nm]');
ylim([-50 350]); xlim([0, 1100])

% Add dummy scatter plots for legend
h1 = scatter(nan, nan, 60, 'MarkerEdgeColor', edgeColor, 'MarkerFaceColor', day1Color);
h2 = scatter(nan, nan, 60, 'MarkerEdgeColor', edgeColor, 'MarkerFaceColor', day2Color);
h3 = scatter(nan, nan, 60, 'MarkerEdgeColor', edgeColor, 'MarkerFaceColor', day3Color);
h4 = scatter(nan, nan, 60, 'MarkerEdgeColor', edgeColor, 'MarkerFaceColor', day4Color);
h5 = scatter(nan, nan, 60, 'MarkerEdgeColor', edgeColor, 'MarkerFaceColor', day5Color);
h6 = scatter(nan, nan, 60, 'MarkerEdgeColor', edgeColor, 'MarkerFaceColor', day6Color);
% Create legend
legend([h1, h2, h3, h4, h5, h6], 'Day 1', 'Day 2', 'Day 3', 'Day 4', 'Day 5', 'Day 6', 'Location', 'best');

% Effective radius between tip and the sample
figure(); hold on
box on; set(gca,'FontSize', 18, 'Linewidth', 1.5);
scatter(CsInden_mean, CsEffectiveRadius_mean, 60, c, "filled");
errorbar(CsInden_mean, CsEffectiveRadius_mean, CsEffectiveRadius_std, 'o', 'Color', c);
xlabel('Indentation depth [nm]');
ylabel('Effective radius [nm]');
xlim([0 200]);ylim([0 600])

figure(); hold on
box on; set(gca,'FontSize', 18, 'Linewidth', 1.5);
scatter(CsEffectiveRadius_mean, CsEModHertz_mean, 60, c, "filled");
errorbar(CsEffectiveRadius_mean, CsEModHertz_mean, CsEModHertz_std, 'o', 'Color', c);
ylabel('Indentation modulus [kPa]');
xlabel('Effective radius [nm]');
xlim([0 600])

% Find indices where any of the arrays have NaN values
nanIndices = isnan(CsInden_mean) | isnan(CsEffectiveRadius_mean);

% Remove NaN values from each array
CsEffectiveRadius_mean_clean = CsEffectiveRadius_mean(~nanIndices);
CsEffectiveRadius_std_clean = CsEffectiveRadius_std(~nanIndices);
CsInden_mean_clean = CsInden_mean(~nanIndices);
CsInden_std_clean = CsInden_std(~nanIndices);
CsEModHertz_mean_clean = CsEModHertz_mean(~nanIndices); 
CsEModHertz_std_clean = CsEModHertz_std(~nanIndices);


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
