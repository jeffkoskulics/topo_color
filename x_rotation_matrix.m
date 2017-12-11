function Rx = x_rotation_matrix(dTheta)
%rotation_matrix creates the following rotation matrix
% Rx = [ 1      0             0       ]
%      [ 0 cos(theta_x) -sin(theta_x) ]
%      [ 0 sin(theta_x) cos(theta_x)  ]

%Reshape dTheta into a list
theta_x = reshape(dTheta,size(dTheta,1)*size(dTheta,2),1);
numPoints = length(theta_x);

%Solve for Rx
Rx(:,1) = [ones(numPoints,1);zeros(3*numPoints,1);cos(theta_x); ...
    sin(theta_x);zeros(numPoints,1);-sin(theta_x);cos(theta_x)];

%Shape into matrix sized n x m x 3 x 3
Rx = reshape(Rx,size(dTheta,1),size(dTheta,2),3,3);
