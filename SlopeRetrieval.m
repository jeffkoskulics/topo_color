function [u_sn] = SlopeRetrieval(data,u_sn_previous,previous_surface, ...
    u_fp,X,Y,PL,spectra_cube)
%% %%% RETRIEVAL STEP 1: SURFACE NORMALS %%%
% The surface normal calculation is, in itself, a two part calculation. The
% x-rotation is applied first, followed by the y-rotation.

%Set constants
theta = .001;

%% X-Rotation
%First apply a very small rotation about the x-axis
u_sn_x = rotate_normals_around_x_axis(u_sn_previous,theta);
    
%Create the model images for the unrotated and x-rotated surface normals
noro_model = forward_model(u_sn_previous,previous_surface,u_fp,1,X,Y,PL,spectra_cube);
xro_model = forward_model(u_sn_x,previous_surface,u_fp,2,X,Y,PL,spectra_cube);
   
%Solve for the rotation angles and create a rotation matrix  
dTheta = solve_dTheta(data,noro_model,xro_model,theta);
Rx = x_rotation_matrix(dTheta);

u_sn_temp = rotate_surface_normals(Rx,u_sn_previous);

%% Y-Rotation    
%Apply a very small rotation about the y-axis
u_sn_y = rotate_normals_around_y_axis(u_sn_temp,theta);

%Create the model images for the unrotated and y-rotated surface normals
noro_model = forward_model(u_sn_temp,previous_surface,u_fp,1,X,Y,PL,spectra_cube);
yro_model = forward_model(u_sn_y,previous_surface,u_fp,2,X,Y,PL,spectra_cube);

%Solve for the rotation angles and create a rotation matrix  
dTheta = solve_dTheta(data,noro_model,yro_model,theta);
Ry = y_rotation_matrix(dTheta);

u_sn = rotate_surface_normals(Ry,u_sn_temp);