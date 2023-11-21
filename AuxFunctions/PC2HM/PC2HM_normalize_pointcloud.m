function [normalized_pc, scale, offset] = PC2HM_normalize_pointcloud(pointcloud, varargin)
    % Get the maximum and minimum coordinates of the point cloud
    x_max = max(pointcloud(:, 1));
    y_max = max(pointcloud(:, 2));
    z_max = max(pointcloud(:, 3));
    x_min = min(pointcloud(:, 1));
    y_min = min(pointcloud(:, 2));
    z_min = min(pointcloud(:, 3));

    % Calculate the scale and offset values for normalization
    scale = [x_max-x_min, y_max-y_min, z_max-z_min];
    offset = [x_min, y_min, z_min];

    % Check if 'KeepAspects' option is passed
    if nargin > 1 && strcmpi(varargin{1},'KeepAspects')
        % Keep aspect ratios the same by dividing the y and z dimensions by the x dimension scale
        Aspectscale(2) = scale(2) / scale(1);
        Aspectscale(3) = scale(3) / scale(1);
        Aspectscale(1) = 1;
        scale = scale./Aspectscale;
    end

    % Normalize the point cloud
    normalized_pc = (pointcloud - offset) ./ scale;
end
