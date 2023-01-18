function PC2HM_visualizePartition(subsets, boundaries, varargin)

% Parse optional Name-Value arguments
p = inputParser;
addParameter(p, 'Scale', []);
addParameter(p, 'Figure', []);
parse(p, varargin{:});
Figure = p.Results.Figure;
% Initialize a figure
gcf = Figure;
hold on;

% Loop through each subset and boundary
for i = 1:length(subsets)
    
    % Plot the points in the subset
    scatter3(subsets{i}(:,1), subsets{i}(:,2), subsets{i}(:,3), '.');
    
    % Plot the boundary cube
    minCoords = boundaries{i}(1,:);
    maxCoords = boundaries{i}(2,:);
    plot3([minCoords(1), maxCoords(1), maxCoords(1), minCoords(1), minCoords(1)], ...
          [minCoords(2), minCoords(2), maxCoords(2), maxCoords(2), minCoords(2)], ...
          [minCoords(3), minCoords(3), minCoords(3), minCoords(3), minCoords(3)], 'k-');
    plot3([minCoords(1), maxCoords(1), maxCoords(1), minCoords(1), minCoords(1)], ...
          [minCoords(2), minCoords(2), maxCoords(2), maxCoords(2), minCoords(2)], ...
          [maxCoords(3), maxCoords(3), maxCoords(3), maxCoords(3), maxCoords(3)], 'k-');
    plot3([minCoords(1), minCoords(1), minCoords(1)], ...
          [minCoords(2), minCoords(2), maxCoords(2)], ...
          [minCoords(3), maxCoords(3), maxCoords(3)], 'k-');
    plot3([maxCoords(1), maxCoords(1), maxCoords(1)], ...
          [minCoords(2), minCoords(2), maxCoords(2)], ...
          [minCoords(3), maxCoords(3), maxCoords(3)], 'k-');
end

% Set the correct physical scale if specified
if ~isempty(p.Results.Scale)
    axis equal;
    axis(p.Results.Scale);
end

end
