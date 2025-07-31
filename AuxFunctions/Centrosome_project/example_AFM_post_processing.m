m = 48; 

ChannelContactHeight = E.FM{m}.get_channel('Contact Height');
figure('name','Contact Height','visible',show_fig); hold on
imagesc(ChannelContactHeight.Image.*1e9); axis image;
c = colorbar; c.Location = 'eastoutside'; c.Label.String = 'Height [nm]';
set(gca,'FontSize', 19, 'Linewidth', 1.5); axis off; colormap("hot")


ChannelContactHeight = E.FM{m}.get_channel('Contact Height Smoothed');
figure('name','Contact Height Smoothed','visible',show_fig); hold on
imagesc(ChannelContactHeight.Image.*1e9); axis image;
c = colorbar; c.Location = 'eastoutside'; c.Label.String = 'Height [nm]';
set(gca,'FontSize', 19, 'Linewidth', 1.5); axis off; colormap("hot")


figure('name', 'Topography angle', 'visible', show_fig);
imagesc(flipud(AngleImage)); axis image; c = colorbar;
c.Location = 'eastoutside'; c.Label.String = 'Angle [rad]'; % 90 deg flat
set(gca,'FontSize', 19, 'Linewidth', 1.5); axis off


 % Angle-based segmentation
 T2 = 1.45; % 83.0788802939694 deg % multithresh(AngleImage);
 AngleCsBW = imbinarize(AngleImage,T2).*erodedBW2; % Exclude flat areas outside centrosome
 figure('name', 'Topography angle segmented','visible',show_fig)
 hold on;
 imagesc(AngleCsBW); axis image;  axis off
 set(gca,'FontSize', 19, 'Linewidth', 1.5);


 Channel = E.FM{m1}.get_channel('Contact Height');
 pxSize1 = Channel.ScanSizeX/Channel.NumPixelsX; % Pixel size in meters
 physWidth1 = size(ChannelContactHeight,2) * pxSize1; % Physical width in meters
 physHeight1 = size(ChannelContactHeight,1) * pxSize1; % Physical height in meters

figure('name', 'Topography angle segmented','visible',show_fig)
hold on;
imagesc(AngleCsBW); 
axis image;  
axis off
set(gca,'FontSize', 19, 'Linewidth', 1.5);

% Calculate scale bar parameters
scalebar_length = 1e-6; % 1 micron in meters
scalebar_pixels = round(scalebar_length / pxSize1); % Convert to pixels

% Position the scale bar in bottom right (10% inset from edges)
x_pos = size(AngleCsBW,2) - scalebar_pixels - round(0.1*size(AngleCsBW,2));
y_pos = size(AngleCsBW,1) - round(0.1*size(AngleCsBW,1));

% Draw the scale bar
rectangle('Position', [x_pos, y_pos, scalebar_pixels, 5], ...
          'FaceColor', 'w', 'EdgeColor', 'w', 'LineWidth', 2);

