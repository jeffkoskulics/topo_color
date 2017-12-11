function dTheta = solve_dTheta(data,noro_model,ro_model,theta)
% Least Squares Solution
% dtheta = (J'J)^-1*J'*dM = J'dotdM/mag(J)
% Expected results:
%   Large changes in the y component of normals will create 
%       a large dTheta_x
%   Large changes in the x component of normals will create
%       a large dTheta_y

%Create the Jacobian matrix
%J(n,m,:,:) = [ dRed/dTheta_x_or_y   ]
%             [ dGreen/dTheta_x_or_y ]
%x_ro is a rotation about the x-axis, so it has a dTheta_x
%y_ro is a rotation about the y-axis, so it has a dTheta_y
J(:,:,:) = ro_model(:,:,1:3) - noro_model(:,:,1:3);
J = J/theta;

deltaColor = data - noro_model;
dTheta = dot(J,deltaColor(:,:,1:3),3)./dot(J,J,3);

%Limit rotations to a maximum angle
%max_angle = 30 * pi / 180;
%dTheta(lt(atan(dTheta),-max_angle)) = -max_angle;
%dTheta(gt(atan(dTheta),max_angle))  =  max_angle;
