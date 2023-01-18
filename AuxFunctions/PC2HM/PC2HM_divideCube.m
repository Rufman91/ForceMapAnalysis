
function [subMinCoords, subMaxCoords] = PC2HM_divideCube(minCoords, maxCoords, factor)

% Calculate the center of the cube
center = (minCoords + maxCoords) / 2;

% Calculate the dimensions of the smaller cubes
subDim(1:2) = (maxCoords(1:2) - minCoords(1:2)) / factor;
subDim(3) = maxCoords(3) - minCoords(3);
% Initialize variables to store the minimum and maximum coordinates of the smaller cubes
subMinCoords = zeros(factor^2, 3);
subMaxCoords = zeros(factor^2, 3);

% Calculate the minimum and maximum coordinates for each of the smaller cubes
idx = 1;
for i = 1:factor
    for j = 1:factor
        subMinCoords(idx,:) = minCoords + [subDim(1) * (i-1), subDim(2) * (j-1), 0];
        subMaxCoords(idx,:) = subMinCoords(idx,:) + subDim;
        idx = idx + 1;
    end
end

end