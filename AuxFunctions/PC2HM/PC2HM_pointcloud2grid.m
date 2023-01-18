function gridArray = PC2HM_pointcloud2grid(points, resolution)
    % Create an empty grid array
    minCoords = min(points);
    maxCoords = max(points);
    gridArray = zeros(resolution, resolution);
    % Fill in the array with the Z values of the points
    for i = 1:length(points)
        x = round(1 + (points(i, 1) - minCoords(1)) * (resolution-1) / (maxCoords(1) - minCoords(1)));
        y = round(1 + (points(i, 2) - minCoords(2)) * (resolution-1) / (maxCoords(2) - minCoords(2)));
        gridArray(x, y) = points(i, 3);
    end
end