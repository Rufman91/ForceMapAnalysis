% This script plots two AFM images with different pixel dimensions 
% at their true physical scale with consistent 1 micron scale bars

m1 = 13; % 44
m2 = 43; % 43 47
show_fig = 'on'; 

% Get first image (68x68)
% Channel1 = E.FM{m1}.get_channel('Contact Height Smoothed');
Channel1 = E.FM{m1}.get_channel('Indentation Modulus Hertz (02)');
% img1 = Channel1.Image.*1e9; % Convert to nm
img1 = Channel1.Image*1e-3; img1(imag(img1) ~= 0) = inf; img1(isnan(img1)) = inf; 
pxSize1 = Channel1.ScanSizeX/Channel1.NumPixelsX; % Pixel size in meters
physWidth1 = size(img1,2) * pxSize1; % Physical width in meters
physHeight1 = size(img1,1) * pxSize1; % Physical height in meters

% Get second image (88x128)
% Channel2 = E.FM{m2}.get_channel('Contact Height Smoothed');
Channel2 = E.FM{m2}.get_channel('Indentation Modulus Hertz (02)');
% img2 = Channel2.Image.*1e9; % Convert to nm
img2 = Channel2.Image*1e-3; img2(imag(img2) ~= 0) = inf; img2(isnan(img2)) = inf; 
pxSize2 = Channel2.ScanSizeX/Channel2.NumPixelsX; % Pixel size in meters
physWidth2 = size(img2,2) * pxSize2; % Physical width in meters
physHeight2 = size(img2,1) * pxSize2; % Physical height in meters

% Create figure
figure('name','Contact Height Comparison','visible',show_fig, ...
       'Units','normalized','Position',[0.1 0.1 0.8 0.6]);

% Plot first image
subplot(1,2,1);
imagesc([0 physWidth1], [0 physHeight1], img1);
axis image; axis off;
% caxis([0 max([img1(:); img2(:)])]); % Shared color scale
caxis([0 350])
% colormap hot; colorbar('northoutside'); 
set(gca, 'FontSize', 19, 'LineWidth', 1.5);
% title(sprintf('Image 1 (%.1f μm × %.1f μm)',physWidth1*1e6,physHeight1*1e6));

% Plot second image with correct aspect ratio
subplot(1,2,2);
imagesc([0 physWidth2], [0 physHeight2], img2);
axis image; axis off;
% caxis([0 max([img1(:); img2(:)])]); % Shared color scale
caxis([0 350])
colormap hot; colorbar('eastoutside');
set(gca, 'FontSize', 19, 'LineWidth', 1.5);
% title(sprintf('Image 2 (%.1f μm × %.1f μm)',physWidth2*1e6,physHeight2*1e6));

% Add consistent scale bar (1 micron)
scaleBarLength = 1e-6; % 1 micron in meters
scaleBarHeight = 0.05 * scaleBarLength; % Height of scale bar

% For first image
subplot(1,2,1);
hold on;
rectangle('Position',[0.8*physWidth1-scaleBarLength, 0.85*physHeight1, ...
                     scaleBarLength, scaleBarHeight], ...
          'FaceColor','k','EdgeColor','k');
% text(0.8*physWidth1-scaleBarLength/2, 0.8*physHeight1, ...
%      '1 μm','Color','k','HorizontalAlignment','center');

% For second image
subplot(1,2,2);
hold on;
rectangle('Position',[0.8*physWidth2-scaleBarLength, 0.85*physHeight2, ...
                     scaleBarLength, scaleBarHeight], ...
          'FaceColor','k','EdgeColor','k');
% text(0.8*physWidth2-scaleBarLength/2, 0.8*physHeight2, ...
%      '1 μm','Color','k','HorizontalAlignment','center');

