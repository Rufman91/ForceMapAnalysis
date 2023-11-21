function [subsets1, boundaries1, subsets2, boundaries2] = PC2HM_partitionPointCloud(points1, values1, points2, values2, MaxSumPerCube, DividingFactor, ExpansionFactor)

% Find the boundary cube around the point cloud
minCoords = min(points1);
maxCoords = max(points1);

% Initialize variables to store the subsets and boundaries
subsets1 = {};
boundaries1 = {};
subsets2 = {};
boundaries2 = {};

% Recursively partition the point cloud into smaller cubes
[subsets1, boundaries1, subsets2, boundaries2] = PC2HM_partitionPointCloudRecursive(...
    points1, values1, points2, values2, minCoords, maxCoords, MaxSumPerCube, DividingFactor, ExpansionFactor,...
    subsets1, boundaries1, subsets2, boundaries2);

end
