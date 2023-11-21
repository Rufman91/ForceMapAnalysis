function [rx, ry, rz] = rotateSphere(x, y, z, center, axis, angle)
    %ROTATESPHERE Rotate a sphere about an arbitrary axis
    % Inputs:
    % x, y, z - Sphere coordinates
    % center - Center of rotation
    % axis - Axis of rotation
    % angle - Angle of rotation (in degrees)

    % Translate the sphere to the origin
    x = x - center(1);
    y = y - center(2);
    z = z - center(3);

    % Normalize the rotation axis
    axis = axis / norm(axis);

    % Convert the angle to radians
    angle = -deg2rad(angle);

    % Set up the rotation matrix
    c = cos(angle);
    s = sin(angle);
    C = 1 - c;
    u = axis(1);
    v = axis(2);
    w = axis(3);
    R = [u^2*C+c, u*v*C-w*s, u*w*C+v*s;
         v*u*C+w*s, v^2*C+c, v*w*C-u*s;
         w*u*C-v*s, w*v*C+u*s, w^2*C+c];

    % Rotate the points
    [m, n] = size(x);
    coords = [x(:), y(:), z(:)] * R';

    % Reshape the rotated coordinates
    rx = reshape(coords(:, 1), m, n) + center(1);
    ry = reshape(coords(:, 2), m, n) + center(2);
    rz = reshape(coords(:, 3), m, n) + center(3);
end
