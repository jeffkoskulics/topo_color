function [u_sn_y] = rotate_normals_around_y_axis(u_sn,theta)
%Rotate surface normals about the y-axis by angle theta

%Rotation matrix
R = [cos(theta)     0           sin(theta)  ;    ...
     0              1           0           ;    ...
     -sin(theta)    0           cos(theta)  ];
 
%Matrix-vector products
u_sn_y = reshape(u_sn,1024*1280,3);
u_sn_y = R*squeeze(u_sn_y(:,:))';
u_sn_y = reshape(u_sn_y',1024,1280,3);