% PLOT_AFM_3DSURFACE Creates publication-quality 3D surface plot of AFM data

%% Extract and prepare data
Channel1 = E.FM{1}.get_channel('Contact Height');
img_um = Channel1.Image * 1e6; % Convert height from meters to nanometers
% Get dimensions from the image itself
[nRows, nCols] = size(img_um);

% Create grid matching the image dimensions
[X,Y] = meshgrid(...
    linspace(0, Channel1.ScanSizeX*1e6, nCols), ...
    linspace(0, Channel1.ScanSizeY*1e6, nRows));

%% Create surface plot with enhanced rendering
surf(X, Y, img_um, 'EdgeColor', 'none', 'FaceColor', 'interp');
axis tight;
daspect([1 1 max(img_um(:))/(mean([Channel1.ScanSizeX Channel1.ScanSizeY])*1e6)]); % Scale z relative to x,y
grid off;

%% Set color mapping
% colormap(parula(256)); % Scientific colormap
% c = colorbar('eastoutside');
% c.Label.String = 'Height [μm]';
% c.Label.FontSize = 17;
% clim = [prctile(img_um(:), 1), prctile(img_um(:), 99)]; % Robust color limits
caxis([0 0.4]);

%% Format axes
% xlabel('[μm]', 'FontSize', 10, 'FontName', 'Helvetica');
% ylabel('[μm]', 'FontSize', 10, 'FontName', 'Helvetica');
% zlabel('[μm]', 'FontSize', 10, 'FontName', 'Helvetica');

%% Lighting for better visualization
material([0.3 0.8 0.2]);
lightangle(-45, 30);
lighting gouraud;
camlight headlight;

set(gca, 'FontSize', 19, 'LineWidth', 1.5);

% view(-30, 60); 
% view(-37.5, 30);  % Classic 3D view that works for most cases
% view(-20, 25);  % More aligned with x-axis
view(-10, 70);  % Higher elevation angle
% view(0, 75);  % Nearly overhead but with height perspective
% view(-60, 40);  % More dramatic perspective
% view(30, 30);  % Different angle but still balanced
rotate3d on;


