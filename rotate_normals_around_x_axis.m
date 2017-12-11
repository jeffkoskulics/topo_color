function [u_sn_x] = rotate_normals_around_x_axis(u_sn,theta)
%Rotate surface normals about the x-axis by angle theta

%Rotation matrix
R = [1  0           0           ;    ...
     0  cos(theta)  -sin(theta) ;    ...
     0  sin(theta)  cos(theta)  ];
 
%Matrix-vector products
u_sn_x = reshape(u_sn,1024*1280,3);
u_sn_x = R*squeeze(u_sn_x(:,:))';
u_sn_x = reshape(u_sn_x',1024,1280,3);