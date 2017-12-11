function Ry = y_rotation_matrix(dTheta)
%rotation_matrix creates the following rotation matrix
% Ry = [ cos(theta_y)   0  sin(theta_y) ]
%      [      0         1       0       ]
%      [ -sin(theta_y)  0  cos(theta_y) ]

%Reshape dTheta into a list
theta_y = reshape(dTheta,size(dTheta,1)*size(dTheta,2),1);
numPoints = length(theta_y);

%Solve for Ry
Ry(:,1) = [cos(theta_y);zeros(numPoints,1);-sin(theta_y); ...
    zeros(numPoints,1);ones(numPoints,1);zeros(numPoints,1); ...
    sin(theta_y);zeros(numPoints,1);cos(theta_y)];

%Shape into matrix sized n x m x 3 x 3
Ry = reshape(Ry,size(dTheta,1),size(dTheta,2),3,3);
