function [output_image] = forward_model(u_sn,r_surf,u_fp,extrapval,...
    X,Y,PL,spectra_cube)
%% Forward_model runs the model on the surface normals.
% An RGB image is output which corresponds to the expected response of the
% system to the given surface normals. Note: If you would like to generate
% surface normals from a r_surf matrix, use the script "surface_normals".
rows = size(u_fp,1);
cols = size(u_fp,2);

%% Tranmitted Ray Unit Vectors
% T = n1/n2*I - (n1/n2*(I*N) + sqrt(1-((n1/n2)^2*(1-(I*N)^2)^2)))*N
%where T is the transmitted vector matrix, I is the incident vector matrix,
%and N is the normal vector matrix
air_index = 1; %n1
water_index = 1.333; %n2h
rel_index = air_index / water_index;

%Calculate constants
fp_dot_sn = dot(u_fp,u_sn,3);
sqrt_consts = (rel_index*fp_dot_sn + sqrt(1-(rel_index^2* ...
    (1-fp_dot_sn.^2).^2)));
%Transmitted unit vectors
u_trans = zeros(rows,cols,3);
for i = 1:3
    u_trans(:,:,i) = rel_index*u_fp(:,:,i) - sqrt_consts.*u_sn(:,:,i);
end
mag = sqrt(dot(u_trans,u_trans,3));
u_trans(:,:,1) = u_trans(:,:,1) ./ mag; %Normalize
u_trans(:,:,2) = u_trans(:,:,2) ./ mag;
u_trans(:,:,3) = u_trans(:,:,3) ./ mag;

%% Fresnel Reflectance
% R = Rs + Rp where the incident angle is defined by the angle between the
% negation of the surface normal and the transmitted ray (assumes both are
% normalized)
theta_i = acos(dot(-u_sn,u_trans,3));

Rs = ((water_index*cos(theta_i) - air_index*sqrt(1 - ((water_index/air_index)*sin(theta_i)).^2))./...
    ((water_index*cos(theta_i) + air_index*sqrt(1 - ((water_index/air_index)*sin(theta_i)).^2)))).^2;

Rp = ((water_index*sqrt(1-((water_index/air_index)*sin(theta_i)).^2) - air_index*cos(theta_i))./...
    (water_index*sqrt(1-((water_index/air_index)*sin(theta_i)).^2) + air_index*cos(theta_i))).^2;

R = (Rs + Rp)/2;

%Transmittance = 1 - Reflectance
T = 1 - R;

%% Bottom Positions
%Bottom intercept solution (tracing the ray from the surface elevation
%point to its intercept on the bottom surface.  Here we're given a
%transmitted vector, and a point on the surface.  We can solve a system of
%parametric equations for z = 0

%x(i,j) = x_surf(i,j) + u_trans(i,j,1)*t
%y(i,j) = y_surf(i,j) + u_trans(i,j,2)*t
%z(i,j) = z_surf(i,j) + u_trans(i,j,3)*t

r_bot = zeros(rows,cols,2);
r_bot(:,:,1) = r_surf(1:rows,1:cols,1) - r_surf(1:rows,1:cols,3).* ...
   (u_trans(:,:,1)./u_trans(:,:,3));
r_bot(:,:,2) = r_surf(1:rows,1:cols,2) - r_surf(1:rows,1:cols,3).* ... 
   (u_trans(:,:,2)./u_trans(:,:,3));

%% Path length through water
%Finds the distance a light ray travels from the point of entrance into the
%water to the point in which it hits the bottom of the tank.
path_length = sqrt((r_surf(1:size(r_surf,1)-1,1:size(r_surf,2)-1,1)-r_bot(:,:,1)).^2 + ...
    (r_surf(1:size(r_surf,1)-1,1:size(r_surf,2)-1,2)-r_bot(:,:,2)).^2 + ...
    (r_surf(1:size(r_surf,1)-1,1:size(r_surf,2)-1,3)).^2);

%% Generate Model Image
%This version uses interp3 to interact with spectrum-created data cube
modelr = interp3(X,Y,PL,spectra_cube(:,:,:,1),r_bot(:,:,1),r_bot(:,:,2),path_length,'*linear',extrapval);
modelg = interp3(X,Y,PL,spectra_cube(:,:,:,2),r_bot(:,:,1),r_bot(:,:,2),path_length,'*linear',extrapval);
modelb = interp3(X,Y,PL,spectra_cube(:,:,:,3),r_bot(:,:,1),r_bot(:,:,2),path_length,'*linear',extrapval);

%Include effects of Fesnel Reflectance
modelr = modelr .* T;
modelg = modelg .* T;
modelb = modelb .* T;

%Find the bad pixels from the previous run and set them to a constant
bad_pixels = lt(u_sn(:,:,3),-1);
modelr(bad_pixels) = extrapval;
modelg(bad_pixels) = extrapval;
modelb(bad_pixels) = extrapval;

output_image(:,:,1) = modelr;
output_image(:,:,2) = modelg;
output_image(:,:,3) = modelb;