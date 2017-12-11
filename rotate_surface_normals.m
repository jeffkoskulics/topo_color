function u_sn = rotate_surface_normals(R,old_u_sn)
%Rotate the surface normals using the rotation matrix around the x or y
%direction.

%We are currently using the new_u_sn = (Rx or Ry)*old_u_sn method.

%Reshape the R martix
D = reshape(R,size(R,1)*size(R,2)*9,1);
numPoints = length(D)/9;

%Make new_u_sn
% form: R*u_sn = [ a b c ] * [ u_x ] = [ ax + by + cz ]
%                [ d e f ]   [ u_y ]   [ dx + ey + fz ]
%                [ h i j ]   [ u_z ]   [ hx + iy + jz ]
U = reshape(old_u_sn,size(old_u_sn,1)*size(old_u_sn,2)*3,1);
P1(:,1) = D(1:numPoints) .* U(1:numPoints) + ...
    D(3*numPoints+1:4*numPoints) .* U(numPoints+1:2*numPoints) + ...
    D(6*numPoints+1:7*numPoints) .* U(2*numPoints+1:3*numPoints);
P2(:,1) = D(numPoints+1:2*numPoints) .* U(1:numPoints) + ...
    D(4*numPoints+1:5*numPoints) .* U(numPoints+1:2*numPoints) + ...
    D(7*numPoints+1:8*numPoints) .* U(2*numPoints+1:3*numPoints);
P3(:,1) = D(2*numPoints+1:3*numPoints) .* U(1:numPoints) + ...
    D(5*numPoints+1:6*numPoints) .* U(numPoints+1:2*numPoints) + ...
    D(8*numPoints+1:9*numPoints) .* U(2*numPoints+1:3*numPoints);
P = [P1;P2;P3];

%Reshape to a n x m x 3 matrix
u_sn = reshape(P,size(old_u_sn,1),size(old_u_sn,2),3);

% SHOULD NOT NEED NORMALIZATION?
%mag = sqrt(dot(u_sn,u_sn,3));
%u_sn(:,:,1) = u_sn(:,:,1) ./ mag; %Normalize
%u_sn(:,:,2) = u_sn(:,:,2) ./ mag;
%u_sn(:,:,3) = u_sn(:,:,3) ./ mag;