function [gridPoints,minCoords, maxCoords] = PC2HM_createGrid(points, resolution)

% Find the boundary cube around the point cloud
minCoords = min(points);
maxCoords = max(points);

% Create a grid of points covering the boundary cube
x = linspace(minCoords(1),maxCoords(1),resolution);
y = linspace(minCoords(2),maxCoords(2),resolution);
[X,Y] = meshgrid(x, y);

% Set the z-value of the grid points to 0
Z = zeros(size(X));

% Combine the x, y, and z values into a single point cloud
gridPoints = [X(:), Y(:), Z(:)];

end
