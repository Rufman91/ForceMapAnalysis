cd(E.ExperimentFolder)
close all

s2 = ' (02)';

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
        % Calculate volume using manual segmentation
        Height{i} = E.FM{i}.get_segment_data_from_channel('Contact Height Smoothed', 'MatchString', 'Seg-02'); % Total centrosome height
        positiveHeight = max(Height{i}, 0); % Treat negative heights as zero
        Volume = sum(positiveHeight) * (E.FM{i}.ScanSizeX/E.FM{i}.NumPixelsX * E.FM{i}.ScanSizeY/E.FM{i}.NumPixelsY); % Total centrosome volume from Seg-02
        Volumes(i) = Volume*1e+18;
        %         % Equivalent radius of a sphere of the same volume
        %         EquivalentRadii_mnl(i) = ((3 * Volumes(i) / (4 * pi))^(1/3))*1000;
    end
end

% Equivalent radius of a sphere of the same volume
EquivalentRadii_mnl = (4.*Volumes./(3*pi)).^(1/3)*1000;

% figure(); hold on
% 
% % Create colormap (one distinct color per ForceMap)
% colors = lines(E.NumForceMaps);
% 
% % Preallocate legend handles
% h_legend = gobjects(E.NumForceMaps, 1);
% 
% % Plot dummy points for legend (invisible)
% for m = 1:E.NumForceMaps
%     h_legend(m) = scatter(nan, nan, 40, ...
%         'MarkerFaceColor', colors(m,:), ...
%         'MarkerEdgeColor', colors(m,:), ...
%         'MarkerFaceAlpha', 0.6, ...
%         'DisplayName', ['ForceMap ' num2str(m)]);
% end
% 
% % Plot each ForceMap's data
% for m = 1:E.NumForceMaps
%     if ismember(m, datavalues)
%     else
%         scatter(all_R2{m}, all_EmodHertz{m}, 40, ...
%             'MarkerFaceColor', colors(m,:), ...
%             'MarkerEdgeColor', colors(m,:), ...
%             'MarkerFaceAlpha', 0.6, ...
%             'HandleVisibility', 'off'); % Exclude from legend
%     end
% end
% 
% % Create legend only using dummy points
% legend(h_legend(~ismember(1:E.NumForceMaps, datavalues)), ...
%     'Location', 'best', 'FontSize', 10, 'NumColumns', 2);
% legend('boxoff');
% 
% xlabel('Hertz Fit R² (02)');
% ylabel('Indentation Modulus (Pa)');
% grid on;
% set(gca, 'FontSize', 12, 'LineWidth', 1.2);
% 
% % Optional: Add R² threshold line if needed
% % R2Thr = 0.99;  % Your threshold value
% xline(R2Thr, '--r', ['R² Threshold = ' num2str(R2Thr)], ...
%     'LineWidth', 1.5, 'LabelOrientation', 'horizontal');
% ylim([0 1.5e7])
% xlim([0.7 1])
% set(gca, 'YScale', 'log')
% 
% 
% figure(); hold on
% 
% % Create colormap (one distinct color per ForceMap)
% colors = lines(E.NumForceMaps);
% 
% % Preallocate legend handles
% h_legend = gobjects(E.NumForceMaps, 1);
% 
% % Plot dummy points for legend (invisible)
% for m = 1:E.NumForceMaps
%     h_legend(m) = scatter(nan, nan, 40, ...
%         'MarkerFaceColor', colors(m,:), ...
%         'MarkerEdgeColor', colors(m,:), ...
%         'MarkerFaceAlpha', 0.6, ...
%         'DisplayName', ['ForceMap ' num2str(m)]);
% end
% 
% % Plot each ForceMap's data
% for m = 1:E.NumForceMaps
%     if ismember(m, datavalues)
%     else
%         scatter(all_PredictiveR2{m}, all_EmodHertz{m}, 40, ...
%             'MarkerFaceColor', colors(m,:), ...
%             'MarkerEdgeColor', colors(m,:), ...
%             'MarkerFaceAlpha', 0.6, ...
%             'HandleVisibility', 'off'); % Exclude from legend
%     end
% end
% 
% % Create legend only using dummy points
% legend(h_legend(~ismember(1:E.NumForceMaps, datavalues)), ...
%     'Location', 'best', 'FontSize', 10, 'NumColumns', 2);
% legend('boxoff');
% 
% xlabel('Hertz Fit R² (02)');
% ylabel('Indentation Modulus (Pa)');
% grid on;
% set(gca, 'FontSize', 12, 'LineWidth', 1.2);
% 
% % Add Predictive R² threshold line
% % PredictiveR2Thr = 0.96;
% xline(PredictiveR2Thr, '--r', ['Predictive R² Threshold = ' num2str(PredictiveR2Thr)], ...
%     'LineWidth', 1.5, 'LabelOrientation', 'horizontal');
% 
% ylim([0 1.5e7]);  % Modulus range
% xlim([0.7 1]);    % Predictive R² range
% set(gca, 'YScale', 'log')
% 
% 
% figure(); hold on
% 
% % Create colormap (one distinct color per ForceMap)
% colors = lines(E.NumForceMaps);
% 
% % Preallocate legend handles
% h_legend = gobjects(E.NumForceMaps, 1);
% 
% % Plot dummy points for legend (invisible)
% for m = 1:E.NumForceMaps
%     h_legend(m) = scatter(nan, nan, 40, ...
%         'MarkerFaceColor', colors(m,:), ...
%         'MarkerEdgeColor', colors(m,:), ...
%         'MarkerFaceAlpha', 0.6, ...
%         'DisplayName', ['ForceMap ' num2str(m)]);
% end
% 
% for m = 1:E.NumForceMaps
%     if ismember(m, datavalues)
%     else
%         scatter(all_R2{m}, all_PredictiveR2{m}, 40, ...
%             'MarkerFaceColor', colors(m,:), ...
%             'MarkerEdgeColor', colors(m,:), ...
%             'MarkerFaceAlpha', 0.6, ...
%             'HandleVisibility', 'off'); % Exclude from legend
%     end
% end
% 
% % Create legend only using dummy points
% legend(h_legend(~ismember(1:E.NumForceMaps, datavalues)), ...
%     'Location', 'best', 'FontSize', 10, 'NumColumns', 2);
% legend('boxoff');
% 
% xlabel('Hertz Fit R² (02)');
% ylabel('Hertz Fit Predictive R² (02)');
% grid on;
% set(gca, 'FontSize', 12, 'LineWidth', 1.2);
% ylim([0.5 1]);
% xlim([0.5 1]);
% 
% % Add Predictive R² threshold line
% yline(PredictiveR2Thr, '--r', ['Predictive R² Threshold = ' num2str(PredictiveR2Thr)], ...
%     'LineWidth', 1.5, 'LabelOrientation', 'horizontal');
% xline(R2Thr, '--r', ['R² Threshold = ' num2str(R2Thr)], ...
%     'LineWidth', 1.5, 'LabelOrientation', 'horizontal');

figure(); 
hold on;

% 1. Set up colormap
num_colors = 256;
cmap = parula(num_colors);
min_radius = min(EquivalentRadii_mnl);
max_radius = max(EquivalentRadii_mnl);

% 2. Create dummy points for legend (one per unique radius or binned radii)
[unique_radii, ~, ic] = unique(EquivalentRadii_mnl);
num_legend_entries = length(EquivalentRadii_mnl);

% Create dummy scatter plots for legend
h_legend = gobjects(num_legend_entries, 1);
for k = 1:num_legend_entries
    if ~isnan(EquivalentRadii_mnl(k))
    norm_radius = (EquivalentRadii_mnl(k) - min_radius)/(max_radius - min_radius);
    color_idx = round(1 + (num_colors-1)*norm_radius);
    
    h_legend(k) = scatter(nan, nan, 40, ...
        'MarkerFaceColor', cmap(color_idx,:), ...
        'MarkerEdgeColor', cmap(color_idx,:), ...
        'MarkerFaceAlpha', 0.6, ...
        'DisplayName', sprintf('%.1f nm', unique_radii(k)));
    end 
end

% 3. Plot actual data (color by radius)
for m = 1:E.NumForceMaps
    if ~ismember(m, datavalues)
        norm_radius = (EquivalentRadii_mnl(m) - min_radius)/(max_radius - min_radius);
        color_idx = round(1 + (num_colors-1)*norm_radius);
        
        scatter(all_R2{m}, all_EmodHertz{m}.*1e-3, 40, ...
            'MarkerFaceColor', cmap(color_idx,:), ...
            'MarkerEdgeColor', cmap(color_idx,:), ...
            'MarkerFaceAlpha', 0.6, ...
            'HandleVisibility', 'off'); % Exclude from legend
    end
end

% 4. Add colorbar and labels
c = colorbar;
c.Label.String = 'Equivalent radius (nm)';
colormap(cmap);
caxis([200 1400]);

xlabel('Hertz Fit R² (02)');
ylabel('Indentation modulus (kPa)');
set(gca, 'FontSize', 12, 'LineWidth', 1.2, 'YScale', 'log');
grid on;
xline(R2Thr, '--r', ['R² Threshold = ' num2str(R2Thr)], ...
    'LineWidth', 1.5, 'LabelOrientation', 'horizontal');

figure(); 
hold on;

% 1. Set up colormap
num_colors = 256;
cmap = parula(num_colors);
min_radius = min(EquivalentRadii_mnl);
max_radius = max(EquivalentRadii_mnl);

% 2. Create dummy points for legend (one per unique radius or binned radii)
[unique_radii, ~, ic] = unique(EquivalentRadii_mnl);
num_legend_entries = length(EquivalentRadii_mnl);

% Create dummy scatter plots for legend
h_legend = gobjects(num_legend_entries, 1);
for k = 1:num_legend_entries
    if ~isnan(EquivalentRadii_mnl(k))
    norm_radius = (EquivalentRadii_mnl(k) - min_radius)/(max_radius - min_radius);
    color_idx = round(1 + (num_colors-1)*norm_radius);
    
    h_legend(k) = scatter(nan, nan, 40, ...
        'MarkerFaceColor', cmap(color_idx,:), ...
        'MarkerEdgeColor', cmap(color_idx,:), ...
        'MarkerFaceAlpha', 0.6, ...
        'DisplayName', sprintf('%.1f nm', unique_radii(k)));
    end 
end

% 3. Plot actual data (color by radius)
for m = 1:E.NumForceMaps
    if ~ismember(m, datavalues)
        norm_radius = (EquivalentRadii_mnl(m) - min_radius)/(max_radius - min_radius);
        color_idx = round(1 + (num_colors-1)*norm_radius);
        
        scatter(all_R2{m}, all_PredictiveR2{m}, 40, ...
            'MarkerFaceColor', cmap(color_idx,:), ...
            'MarkerEdgeColor', cmap(color_idx,:), ...
            'MarkerFaceAlpha', 0.6, ...
            'HandleVisibility', 'off'); % Exclude from legend
    end
end

% 4. Add colorbar and labels
c = colorbar;
c.Label.String = 'Equivalent radius (nm)';
colormap(cmap);
caxis([200 1400]);

xlabel('Hertz fit R²')
ylabel('Hertz fit Predictive R²');
grid off;
ylim([0.8 1]);
xlim([0.8 1]);

% Add Predictive R² threshold line
yline(PredictiveR2Thr, '--r', ['Predictive R² Threshold = ' num2str(PredictiveR2Thr)], ...
    'LineWidth', 1.5, 'LabelOrientation', 'horizontal');
xline(R2Thr, '--r', ['R² Threshold = ' num2str(R2Thr)], ...
    'LineWidth', 1.5, 'LabelOrientation', 'horizontal');




figure('Name', 'Centrosome volume dependence'); 
hold on;
box on; 
set(gca,'FontSize', 19, 'Linewidth', 1.5);

% Define colors
orange = [0.8500 0.3250 0.0980]; % MATLAB default orange
blue = [0 0.4470 0.7410]; % MATLAB default blue

% Initialize handles for legend
h_blue = [];
h_orange = [];

% Plot individual data points
for m = 1:length(EquivalentRadii_mnl)
    if ~ismember(m, datavalues)
        % Skip Day 2 data (indices 14-32)
        if m >= 14 && m <= 32
            continue;
        end

        % Get current data and filter out zeros
        current_data = all_EmodHertz{m}.*1e-3;
        non_zero_indices = current_data ~= 0;
        filtered_data = current_data(non_zero_indices);

        % Skip if all data points are zero
        if isempty(filtered_data)
            continue;
        end

        % Determine color based on radius
        if EquivalentRadii_mnl(m) > 500
            point_color = orange;
            h = scatter(EquivalentRadii_mnl(m)*ones(size(filtered_data)), ...
                filtered_data, 40, ...
                'MarkerFaceColor', orange, ...
                'MarkerEdgeColor', orange, ...
                'MarkerFaceAlpha', 0.6);
            h_orange = h; % Store handle for legend
        else
            point_color = blue;
            h = scatter(EquivalentRadii_mnl(m)*ones(size(filtered_data)), ...
                filtered_data, 40, ...
                'MarkerFaceColor', blue, ...
                'MarkerEdgeColor', blue, ...
                'MarkerFaceAlpha', 0.6);
            h_blue = h; % Store handle for legend
        end
    end
end

% Add labels and limits
ylabel('Indentation modulus [kPa]');
xlabel('Centrosome equivalent radius [nm]');
xlim([0 1500]); 

% Create legend using representative points
if ~isempty(h_blue) && ~isempty(h_orange)
    legend([h_blue(1), h_orange(1)], {'Radius ≤ 500 nm', 'Radius > 500 nm'}, ...
           'Location', 'best');
elseif ~isempty(h_blue)
    legend(h_blue(1), {'Radius ≤ 500 nm'}, 'Location', 'best');
elseif ~isempty(h_orange)
    legend(h_orange(1), {'Radius > 500 nm'}, 'Location', 'best');
end
legend box off;

figure(); 
hold on;

c = [0 0.4470 0.7410]; 
for m = 1:E.NumForceMaps
    if ~ismember(m, datavalues)
        scatter(all_R2{m}, all_PredictiveR2{m}, 40, ...
            'MarkerFaceColor', c, ...
            'MarkerEdgeColor', c, ...
            'MarkerFaceAlpha', 0.6, ...
            'HandleVisibility', 'off'); % Exclude from legend
    end
end

xlabel('Hertz fit R²')
ylabel('Hertz fit Predictive R²');
grid off;
ylim([0.8 1]);
xlim([0.8 1]);

% Add Predictive R² threshold line
yline(PredictiveR2Thr, '--r', ['Predictive R² Threshold = ' num2str(PredictiveR2Thr)], ...
    'LineWidth', 1.5, 'LabelOrientation', 'horizontal');
xline(R2Thr, '--r', ['R² Threshold = ' num2str(R2Thr)], ...
    'LineWidth', 1.5, 'LabelOrientation', 'horizontal');

total_nonzero_pixels = 0;
excluded_nonzero_pixels = 0;

for m = 1:E.NumForceMaps
    if ~ismember(m, datavalues)  % Skip manually excluded maps
        R2_map = all_R2{m};
        PredR2_map = all_PredictiveR2{m};
        
        % Create a mask for non-zero pixels (where at least R2 or PredR2 is non-zero)
        non_zero_mask = (R2_map ~= 0) | (PredR2_map ~= 0);
        
        % Count only non-zero pixels
        total_nonzero_pixels = total_nonzero_pixels + sum(non_zero_mask(:));
        
        % Find excluded pixels (only where non-zero)
        excluded = non_zero_mask & ((R2_map < R2Thr) | (PredR2_map < PredictiveR2Thr));
        excluded_nonzero_pixels = excluded_nonzero_pixels + sum(excluded(:));
    end
end

if total_nonzero_pixels > 0
    exclusion_percentage = (excluded_nonzero_pixels / total_nzero_pixels) * 100;
    fprintf('Exclusion percentage (non-zero pixels only): %.2f%%\n', exclusion_percentage);
else
    fprintf('No non-zero pixels to evaluate.\n');
end