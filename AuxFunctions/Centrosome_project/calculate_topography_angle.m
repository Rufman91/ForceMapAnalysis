function [PlaneAngle] = calculate_topography_angle(arr)

% generate test fake data
% create a cloud of random points
%         sz = 4;
%         nPoints = 300;
%         arr = rand(nPoints, 3)*sz;
%         pc = pointCloud(arr)
%         pcshow(pc)

% read data from excel file (Desktop)
% arr = xlsread('Step_scan01_ex.xls');

x=arr(:,1);
y=arr(:,2);
z=arr(:,3);
DM = [x, y, ones(size(z))];                             % design Matrix
B = DM\z;                                               % estimate Parameters
[X,Y] = meshgrid(linspace(min(x),max(x),50), linspace(min(y),max(y),50));
Z = B(1)*X + B(2)*Y + B(3)*ones(size(X)); % plane equation

% figure()
% plot3(arr(:,1),arr(:,2),arr(:,3),'.')
% hold on
% meshc(X, Y, Z)
% hold off
% grid on
% xlabel('x(mm)'); ylabel('y(mm)'); zlabel('z(mm)');
% title('Masked plot');
% grid on
% text(-20, 50, 450, sprintf('Z = %.3f\\cdotX %+.3f\\cdotY %+3.0f', B))

% calculate plane angle with the horizontal
v = [B(1), B(2), -1]; % normal vector to plane
v1 = v/norm(v); % unit length vector 
% vh = [0 0 1]; % normal vector to horizontal plane ?
vh = [B(1), B(2), 0]; 
v2 = vh/norm(vh); % unit length vector 
PlaneAngle = atan2(norm(cross(v1,v2)),dot(v1,v2));
% D = rad2deg( PlaneAngle )

end
