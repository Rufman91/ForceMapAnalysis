% clear
% E = Experiment.load;
cd(E.ExperimentFolder)
% close all

msg1 = "Do you need to exclude force maps?";
opts1 = ["Yes" "No"];
choice1 = menu(msg1,opts1);
choice3 = 2; 
s2 = ' ('+opts3(choice3)+')';

% msg2 = "Do you want to apply a subsequent ForceMapAnalysisOptions?";
% opts2 = ["Yes" "No"];
% choice2 = menu(msg2,opts2);

% if choice2 == 1
%     msg3 = "Which ForceMapAnalysisOptions do you want to apply?";
%     opts3 = ["01" "02" "03" "04" "05"];
%     choice3 = menu(msg3,opts3);
%     s2 = ' ('+opts3(choice3)+')';
% else
%     s2 = '';
%     choice3 = []; 
% end

fileID = fopen('Deconvolution_include.txt', 'r'); % Open the file

deconvAllValues = []; % Variable to store all values (including those with '?')
deconvIntegerValues = []; % Variable to store only integer values

% Read the file line by line
while ~feof(fileID)
    line = fgetl(fileID); % Read each line as a string
    
    % Extract the numeric part from the string
    numericValue = sscanf(line, '%d');
    
    % Check if the line contains a number and is an integer
    if ~isempty(numericValue)
        deconvAllValues = [deconvAllValues; numericValue];  % Add to the list of all values
        
        % If the line does not contain special characters (like '?'), add to integer values
        if isempty(regexp(line, '\?', 'once')) && isempty(regexp(line, 'Half previous height', 'once'))
            deconvIntegerValues = [deconvIntegerValues; numericValue];
        end
    end
end

fclose(fileID); % Close the file

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

% Refine the variables by removing excluded values
refinedDeconvAllValues = setdiff(deconvAllValues, datavalues); % Exclude from all values
refinedDeconvIntegerValues = setdiff(deconvIntegerValues, datavalues); % Exclude from integer values

clear Height
clear Volumes
clear EquivalentRadii
clear Height_deconv
clear Volumes_deconv
clear EquivalentRadii_deconv
clear CsEModHertz_mean
clear CsEModHertz_std

for i = 1:E.NumForceMaps
    if ~ismember(i, refinedDeconvIntegerValues)
        % Skip evaluation round
    else
        cd(E.ForceMapFolders{i,1})
        load(strcat('Processed',s2,'.mat'))
        CsEModHertz_data{i} = CsEModHertz(:);
        CsEModHertz_mean(i) = mean(CsEModHertz(:), 'omitnan').*1e-3;
        CsEModHertz_std(i) = std(CsEModHertz(:), 'omitnan').*1e-3;
        CsFlatHeight_mean(i) = mean(CsFlatHeight(:),'omitnan').*1e9;
        CsFlatPrctile_data(i) = CsFlatPrctile*1e9;
        CsInden_mean(i) = mean(CsFlatInden(:),'omitnan').*1e9;

        % Calculate volume using manual segmentation
        Height{i} = E.FM{i}.get_segment_data_from_channel('Contact Height Smoothed', 'MatchString', 'Seg-02'); % Total centrosome height
        positiveHeight = max(Height{i}, 0); % Treat negative heights as zero
        Volume = sum(positiveHeight) * (E.FM{i}.ScanSizeX/E.FM{i}.NumPixelsX * E.FM{i}.ScanSizeY/E.FM{i}.NumPixelsY); % Total centrosome volume from Seg-02
        Volumes(i) = Volume*1e+18;

        % Calculate volume from deconvoluted channel
        Height_deconv{i} = E.FM{i}.get_segment_data_from_channel('Deconvoluted', 'MatchString', 'Seg-02'); % Total centrosome height
        positiveHeight_deconv = max(Height_deconv{i}, 0); % Treat negative heights as zero
        Volume_deconv = sum(positiveHeight_deconv) * (E.FM{i}.ScanSizeX/E.FM{i}.NumPixelsX * E.FM{i}.ScanSizeY/E.FM{i}.NumPixelsY); % Total centrosome volume from Seg-02
        Volumes_deconv(i) = Volume_deconv*1e+18;
    end
end

CsEModHertz_mean(CsEModHertz_mean == 0) = NaN;
CsEModHertz_std(CsEModHertz_std == 0) = NaN;
CsFlatHeight_mean(CsFlatHeight_mean == 0) = NaN;
CsFlatPrctile_data(CsFlatPrctile_data == 0) = NaN; 
CsInden_mean(CsInden_mean == 0) = NaN; 
Volumes(Volumes == 0) = NaN; 
Volumes_deconv(Volumes_deconv == 0) = NaN; 

% Equivalent radius of a sphere of the same volume
EquivalentRadii = (4.*Volumes./(3*pi)).^(1/3)*1000; 
EquivalentRadii_deconv = (4.*Volumes_deconv./(3*pi)).^(1/3)*1000; 

c = [166/255 54/255 3/255];
% c = [230/255 85/255 13/255];
% c = [253/255 190/255 133/255]; 

figure('name', 'Centrosome volume deconvolution'); hold on
box on; set(gca,'FontSize', 18, 'Linewidth', 1.5);
scatter(EquivalentRadii, EquivalentRadii_deconv, 60, c, "filled");
ylabel('Deconvoluted equivalent radius [nm]'); 
xlabel('Equivalent radius [nm]');

% Linear regression fit
validIdx = ~isnan(EquivalentRadii) & ~isnan(EquivalentRadii_deconv); % Remove NaN values from both variables
EquivalentRadii_valid = EquivalentRadii(validIdx);
EquivalentRadii_deconv_valid = EquivalentRadii_deconv(validIdx);
[p, S] = polyfit(EquivalentRadii_valid, EquivalentRadii_deconv_valid, 1);  % Linear fit
yfit = polyval(p, EquivalentRadii_valid);  % Predicted values based on fit
correlation_coefficient = corr(EquivalentRadii_valid', EquivalentRadii_deconv_valid', 'Type', 'Pearson');

% Calculate the confidence intervals for the slope and intercept
[fit_ci, ~] = polyconf(p, EquivalentRadii_valid, S, 'alpha', 0.05); % 95% CI for the fit

% Calculate R-squared
SS_res = sum((EquivalentRadii_deconv_valid - yfit).^2);  % Sum of squares of residuals
SS_tot = sum((EquivalentRadii_deconv_valid - mean(EquivalentRadii_deconv_valid)).^2);  % Total sum of squares
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


figure('name', 'Centrosome volume dependence'); hold on
box on; set(gca,'FontSize', 18, 'Linewidth', 1.5);
scatter(EquivalentRadii_deconv, CsEModHertz_mean, 60, c, "filled");
errorbar(EquivalentRadii_deconv, CsEModHertz_mean, CsEModHertz_std, 'o', 'Color', c);
ylabel('Indentation modulus [kPa]');
xlabel('Centrosome equivalent radius [nm]');
xlim([0 1600]); ylim([-50 350])

RobustFit = false;

% Find indices where any of the arrays have NaN values
nanIndices = isnan(EquivalentRadii_deconv) | isnan(CsEModHertz_mean);

% Remove NaN values from each array
EquivalentRadii_clean = EquivalentRadii_deconv(~nanIndices);
CsEModHertz_mean_clean = CsEModHertz_mean(~nanIndices);
addCorrelationInfo(EquivalentRadii_clean, CsEModHertz_mean_clean, RobustFit);

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

