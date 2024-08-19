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
        %         CsInden_mean(i) = mean(CsFlatInden(:),'omitnan').*1e9;
        %         CsRadiusXY_data(i) = CsRadiusXY;
        %         CsAspectRatio(i) = mean(CsFlatHeight(:),'omitnan')/(CsRadiusXY*2);
        %         CsFlatArea_data(i) = CsFlatArea;
        CsVolume_data(i) = CsVolume*1e+18; % From automatic segmentation
        CsVolumeSphereCap_data(i) = CsVolumeSphereCap*1e+18;

        Height{i} = E.FM{i}.get_segment_data_from_channel('Contact Height Smoothed', 'MatchString', 'Seg-02'); % Total centrosome height

        % Calculate volume (treat negative heights as zero)
        positiveHeight = max(Height{i}, 0);
        Volume = sum(positiveHeight) * (E.FM{i}.ScanSizeX/E.FM{i}.NumPixelsX * E.FM{i}.ScanSizeY/E.FM{i}.NumPixelsY); % Total centrosome volume from Seg-02
        Volumes(i) = Volume*1e+18;

        % Equivalent radius of a sphere of the same volume
        EquivalentRadii(i) = ((3 * Volumes(i) / (4 * pi))^(1/3))*1000;
    end
end

CsEModHertz_mean(CsEModHertz_mean == 0) = NaN;
CsEModHertz_std(CsEModHertz_std == 0) = NaN;
CsFlatHeight_mean(CsFlatHeight_mean == 0) = NaN;
CsFlatPrctile_data(CsFlatPrctile_data == 0) = NaN; 
CsVolume_data(CsVolume_data == 0) = NaN;
CsVolumeSphereCap_data(CsVolumeSphereCap_data == 0) = NaN;
Volumes(Volumes == 0) = NaN; 
% Equivalent radius of a sphere of the same volume
EquivalentRadii = (4.*Volumes./(3*pi)).^(1/3)*1000; % Use here Volumes (manual) or CsVolume_data (auto)

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
    c =  [152/255 78/255 163/255]; % Thin film
elseif choice3 == 2
    c = [255/255 127/255 0/255]; % Topography
    elseif choice3 == 3
    c = [77/255 175/255 74/255]; % Thin film + topography
end

figure('name', 'Volume dependence'); hold on
box on; set(gca,'FontSize', 18, 'Linewidth', 1.5);
for i = 1:E.NumForceMaps
    scatter(EquivalentRadii(i), CsEModHertz_mean(i), 50, c, "filled");
end
ylabel('Indentation modulus [kPa]');
xlabel('Centrosome equivalent radius [nm]');
xlim([0 1600]); ylim([-50 350])

figure('name', 'Mean height dependence'); hold on
box on; set(gca,'FontSize', 18, 'Linewidth', 1.5);
for i = 1:E.NumForceMaps
    scatter(CsFlatHeight_mean(i), CsEModHertz_mean(i), 50, c, "filled");
end
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
% 
% % Plot fit with data.
% h = plot( fitresult, xData, yData );
ylabel('Indentation modulus [kPa]');
xlabel('Centrosome mean height [nm]');
ylim([-50 350]); xlim([0, 1100])
% legend boxoff

figure('name', 'Maximum height vs. volume'); hold on
box on; set(gca,'FontSize', 18, 'Linewidth', 1.5);
for i = 1:E.NumForceMaps
    scatter(EquivalentRadii(i), CsFlatPrctile_data(i), 50, [0 0.4470 0.7410], "filled");
end
ylabel('Centrosome maximum height [nm]'); 
xlabel('Centrosome equivalent radius [nm]');
ylim([0, 1100]); xlim([0 1800])
 
% figure('name', 'Compression vs. maximum height'); hold on
% box on; set(gca,'FontSize', 16, 'Linewidth', 1.5);
% for i = 1:E.NumForceMaps
%         scatter(CsFlatPrctile_data(i), (CsInden_mean(i)/CsFlatHeight_mean(i))*100, 50, [0.4940 0.1840 0.5560], "filled");  
% end
% xlabel('Centrosome height [nm]');
% ylabel('Compression [%]')


RobustFit = false;

% Find indices where any of the arrays have NaN values
nanIndices = isnan(EquivalentRadii) | isnan(CsEModHertz_mean) | isnan(CsEModHertz_std);

% Remove NaN values from each array
EquivalentRadii_clean = EquivalentRadii(~nanIndices);
CsEModHertz_mean_clean = CsEModHertz_mean(~nanIndices);
CsEModHertz_std_clean = CsEModHertz_std(~nanIndices);

RobustFit = false;

figure;
errorbar(EquivalentRadii_clean, CsEModHertz_mean_clean, CsEModHertz_std_clean, 'o');
xlabel('Centrosome equivalent radius');
ylabel('Mean EMod K1000+Thin Film'); % Simple
title('Mean EMod K1000+Thin Film vs. Centrosome volume');
addCorrelationInfo(EquivalentRadii_clean, CsEModHertz_mean_clean, RobustFit);
box on; set(gca,'FontSize', 18, 'Linewidth', 1.5);
ylim([-50 350]); xlim([0 1600])

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

