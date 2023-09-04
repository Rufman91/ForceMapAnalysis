function [XT, YT, ZT] = transform_3D_PC_coordinates(X, Y, Z, ShiftX, ShiftY, Angle)
    
    % Calculate the center point of the point cloud
    centerX = min(X) + range(X)/2;
    centerY = min(Y) + range(Y)/2;
    
    % Translate point cloud to origin centered on the mean
    X_centered = X - centerX;
    Y_centered = Y - centerY;
    Z_centered = Z; % Z remains unchanged
    
    % Convert angle to radians
    theta_rad = deg2rad(Angle);
    
    % Rotation matrix for rotation around Z-axis
    Rz = [cos(theta_rad), -sin(theta_rad), 0;
          sin(theta_rad),  cos(theta_rad), 0;
          0,               0,              1];
    
    % Form the translated point cloud matrix (3 x N, where N is the number of points)
    pointCloudMatrix = [X_centered'; Y_centered'; Z_centered'];

    % Apply the rotation to the entire point cloud matrix
    rotatedPointCloud = Rz * pointCloudMatrix;
    
    % Translate the point cloud back
    X_rotated_translated = rotatedPointCloud(1, :) + centerX;
    Y_rotated_translated = rotatedPointCloud(2, :) + centerY;
    Z_rotated_translated = rotatedPointCloud(3, :);  % Z remains unchanged
    
    % Shift the point cloud
    XT = X_rotated_translated' + ShiftX;
    YT = Y_rotated_translated' + ShiftY;
    ZT = Z_rotated_translated';  % Z remains unchanged
end
