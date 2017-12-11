function [u_fp, r_surf, X, Y, flat_fitting ] = flat_fit(flat,z_fp,fp_res...
    ,nadir,r_ap)
%% Script to prepare the flat surface data to be used in running
% interp2 during the forward model calculations. A model is used to
% fit the raw data, and the matrices are output in a form used by 
% interp2: color_matrix(i,j,:) = (x_pos,y_pos,color_value);

%% Focal plane position for each pixel
%r_fp(i,j,:) = [x y z], is calculated using the following logic
%x = (j-nadir(2))*fp_res , y = (i-nadir(1))*fp_res , z = z_fp
rows = size(flat,1);
cols = size(flat,2);
r_fp = zeros(rows,cols,3);
r_fp(:,:,1) = repmat((((1:1:cols)-nadir(2))*fp_res),rows,1);
r_fp(:,:,2) = repmat((((1:1:rows)-nadir(1))*fp_res)',1,cols);
r_fp(:,:,3) = repmat(z_fp,rows,cols);

%% Incident Vector 
%located on focal plane pointing along the line from the 
%camera aperature to the pixel location.
%u_fp(i,j,:) = fp_positions(i,j,:) - [0,0,z_aperture];
u_fp = zeros(rows,cols,3);
u_fp(:,:,1) = r_fp(:,:,1) - r_ap(1);
u_fp(:,:,2) = r_fp(:,:,2) - r_ap(2);
u_fp(:,:,3) = r_fp(:,:,3) - r_ap(3);
mag = sqrt(dot(u_fp,u_fp,3));
u_fp(:,:,1) = u_fp(:,:,1) ./ mag; %Normalize
u_fp(:,:,2) = u_fp(:,:,2) ./ mag;
u_fp(:,:,3) = u_fp(:,:,3) ./ mag;


%% Flat surface at focal plane elevation
r_surf = zeros(rows+1,cols+1,3);
r_surf(:,:,1) = repmat((((1:1:cols+1)-nadir(2))*fp_res),rows+1,1);
r_surf(:,:,2) = repmat((((1:1:rows+1)-nadir(1))*fp_res)',1,cols+1);
r_surf(:,:,3) = repmat(z_fp,rows+1,cols+1);

%Define new size
rows = size(u_fp,1);
cols = size(u_fp,2);

u_sn = -1*cross((r_surf(2:rows+1,1:cols,:)-r_surf(1:rows,1:cols,:)),...
    (r_surf(1:rows,2:cols+1,:)-r_surf(1:rows,1:cols,:)),3);
mag = sqrt(dot(u_sn,u_sn,3));
u_sn(:,:,1) = u_sn(:,:,1) ./ mag; %Normalize
u_sn(:,:,2) = u_sn(:,:,2) ./ mag;
u_sn(:,:,3) = u_sn(:,:,3) ./ mag;

%% Tranmitted Ray Unit Vectors
% T = n1/n2*I - (n1/n2*(I*N) + sqrt(1-((n1/n2)^2*(1-(I*N)^2)^2)))*N
%where T is the transmitted vector matrix, I is the incident vector matrix,
%and N is the normal vector matrix
air_index = 1; %n1
water_index = 1.334; %n2h
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
%% Get flat surface fitting
% A quadratic polynomial model is used to fit the color values of the 
% bottom surface. This model ensures smoothness of the data.

red = flat(2:2:rows,1:2:cols);
red = reshape(red,numel(red),1);
red_x = r_bot(2:2:rows,1:2:cols,1);
red_x = reshape(red_x,numel(red_x),1);
red_y = r_bot(2:2:rows,1:2:cols,2);
red_y = reshape(red_y,numel(red_x),1);

red_flat = polyfitn([red_x red_y],red,2);

green = flat(1:2:rows,1:2:cols);
green = reshape(green,numel(green),1);
green_x = r_bot(1:2:rows,1:2:cols,1);
green_x = reshape(green_x,numel(green_x),1);
green_y = r_bot(1:2:rows,1:2:cols,2);
green_y = reshape(green_y,numel(green_y),1);

green_flat = polyfitn([green_x green_y],green,2);

blue = flat(1:2:rows,2:2:cols);
blue = reshape(blue,numel(blue),1);
blue_x = r_bot(1:2:rows,2:2:cols,1);
blue_x = reshape(blue_x,numel(blue_x),1);
blue_y = r_bot(1:2:rows,2:2:cols,2);
blue_y = reshape(blue_y,numel(blue_y),1);

blue_flat = polyfitn([blue_x blue_y],blue,2);

%% Define the plaid grid to be used for interp2

%Define the (approx) number of equally spaced points you want to use
% in interp2. The nearest whole number spacing will be used.

%The physical bounds of x and y color filter data are:
% x: -.2336 to .2340
% y: -.1869 to .1872
x_start = -.5;
x_end = .5;
x_numpoints = 100;

y_start = -.5;
y_end = .5;
y_numpoints = 100;

X = x_start:(x_end-x_start)/x_numpoints:x_end;
Y = y_start:(y_end-y_start)/y_numpoints:y_end;

poly_val_locations(:,:,1) = repmat(X,length(Y),1);
poly_val_locations(:,:,2) = repmat(Y',1,length(X));

%% Solve the color model for points in plaid grid
model(:,1)  = polyvaln(red_flat,reshape(poly_val_locations,length(X)*length(Y),2));
model(:,2)  = polyvaln(green_flat,reshape(poly_val_locations,length(X)*length(Y),2));
model(:,3) = polyvaln(blue_flat,reshape(poly_val_locations,length(X)*length(Y),2));
flat_fitting = reshape(model,length(Y),length(X),3);

