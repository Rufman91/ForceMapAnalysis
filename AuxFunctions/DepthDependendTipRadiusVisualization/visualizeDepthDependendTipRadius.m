function visualizeDepthDependendTipRadius(Channel, radiusVec, heightVec)
    % Validate inputs
    if nargin < 3
        error('Three input arguments required: Channel, radiusVec, heightVec');
    end

    % Extracting the data from Channel
    heightMap = Channel.Image - min(Channel.Image,[],'all'); % Using 'Image' as the height map field
    NumPixelsX = Channel.NumPixelsX;
    NumPixelsY = Channel.NumPixelsY;
    ScanSizeX = Channel.ScanSizeX; % in meters
    ScanSizeY = Channel.ScanSizeY; % in meters

    % Creating real-world coordinate grids
    x = linspace(0, ScanSizeX, NumPixelsX);
    y = linspace(0, ScanSizeY, NumPixelsY);
    [XGrid, YGrid] = meshgrid(x, y);

    % Preparing the figure and slider
    f = figure;
    ax = axes('Parent', f);
    slider = uicontrol('Parent', f, 'Style', 'slider', 'Position', [100, 10, 300, 20],...
        'min', 1, 'max', length(radiusVec), 'Value', 1, 'SliderStep', [1/(length(radiusVec)-1) , 1/(length(radiusVec)-1)]);
    fitRangeText = uicontrol('Style', 'text', 'Position', [410, 10, 100, 20], 'String', 'Fit Range: ');
    addlistener(slider, 'ContinuousValueChange', @(hObject, event) updatePlot(hObject, ax, XGrid, YGrid, heightMap, radiusVec, heightVec, fitRangeText));

    % Initial plot
    updatePlot(slider, ax, XGrid, YGrid, heightMap, radiusVec, heightVec, fitRangeText);
end

function updatePlot(hObject, ax, XGrid, YGrid, heightMap, radiusVec, heightVec, fitRangeText)
    % Preserve the current view
    [az, el] = view(ax);

    % Get the current slider value
    idx = round(get(hObject, 'Value'));

    % Update fit range text
    set(fitRangeText, 'String', sprintf('Fit Range: %.2fm', heightVec(idx)));

    % Generate paraboloid on the same grid
    Z = createParaboloid(XGrid, YGrid, radiusVec(idx), max(heightMap(:)));

    % Adjusting the fitted range map
    fittedRangeMap = NaN(size(heightMap));
    tipApexHeight = max(heightMap(:));
    fittedRangeMap(heightMap >= tipApexHeight - heightVec(idx)) = heightMap(heightMap >= tipApexHeight - heightVec(idx));

    % Clear current axes and replot
    cla(ax);
    hold(ax, 'on');
    surf(ax, XGrid, YGrid, heightMap, 'FaceAlpha', 0.3); % AFM tip with some transparency
    surf(ax, XGrid, YGrid, Z, 'FaceAlpha', 0.5); % Paraboloid with some transparency
    surf(ax, XGrid, YGrid, fittedRangeMap, 'FaceColor', 'r', 'FaceAlpha', 0.5); % Fitted range indicator

    % Adjust visualization properties
    hold(ax, 'off');
    view(ax, [az, el]); % Restore the view
    axis(ax, 'equal'); % Set axes to have equal proportions
    axis(ax, 'tight'); % Fit axes to the data
    xlabel(ax, 'X (m)');
    ylabel(ax, 'Y (m)');
    zlabel(ax, 'Height (m)');
end


function Z = createParaboloid(XGrid, YGrid, radius, maxHeight)
    % Center of the grid
    centerX = mean(XGrid(1,:));
    centerY = mean(YGrid(:,1));

    % Shifting XGrid and YGrid to be centered for paraboloid calculation
    X = XGrid - centerX;
    Y = YGrid - centerY;

    % Paraboloid Equation Z = maxHeight - (x^2 + y^2)/(2*radius)
    Z = maxHeight - ((X.^2 + Y.^2) / (2 * radius));

    % Truncate paraboloid below zero
    Z(Z < 0) = 0;
end

