function elevation = prefactored_solve_global_barge(u_sn)

%%prefactored_solve computes the elevation based on a prefactored Q-less QR
%%factorization.  Inputs needed are weighting matrix W, triangular matrix
%%R, and the measurement vector u_sn

%Identify uncorrupted pixels 
top    = 99; %top row
bottom = 557; %bottom row
left   = 105; %left column
right  = 966; %right column

%Find indices of uncorrupted and corrupted pixels
good_rows = top:bottom;
good_cols = left:right;
good_pixels = zeros(1024,1280);
good_pixels(good_rows,good_cols) = 1;
good_pixels = gt(good_pixels,0.9);
bad_pixels = ~good_pixels;

global W R R_transp A_perm_transp e

tic;
[rows,cols,~] = size(u_sn);

normals(:,1) = reshape(u_sn(:,:,1),rows*cols,1);
normals(:,2) = reshape(u_sn(:,:,2),rows*cols,1);
normals(:,3) = reshape(u_sn(:,:,3),rows*cols,1);

%Set corrupted surface normal to zenith
normals(bad_pixels,1:2) = 0;
normals(bad_pixels,3)   = 1;

%Identify "really" bad slope (those with 60 deg slopes) measurements and 
%set to zero slope
really_bad_pixels = lt(normals(:,3),0.5);
normals(really_bad_pixels,1:2) = 0;  %Assume the slope of bad pixels equals zero
normals(really_bad_pixels,3)   = 1;  %Zero slope has normal pointing straight up

%Organize measurement vector components
b_x = normals(:,1)./normals(:,3);
b_y = normals(:,2)./normals(:,3);
curv_x = zeros((rows+1)*(cols-1),1); 
curv_y = zeros((rows-1)*(cols+1),1); 
z_mean = 0.63*ones((rows+1)*(cols+1),1);

%Determine mean slopes of uncorrupted pixels
mean_slope_x = mean(b_x(good_pixels));
mean_slope_y = mean(b_y(good_pixels));
%Subtract mean slope of uncorrupted pixels from each uncorrupted pixel
b_x(good_pixels) = b_x(good_pixels) - mean_slope_x; 
b_y(good_pixels) = b_y(good_pixels) - mean_slope_y; 

%Assemble measurement vector
b = [b_x; b_y; curv_x; curv_y; z_mean];

tic;
%x = (W*A)\(W*b);

b = W*b;
b = A_perm_transp*b;
y = R_transp\b;
x(e,1) = R\y;

%x(e,1) = R\(R'\(A(:,e)'*W'*W*b));
%r = W*b - W*A*x;
%err(e,1) = R\(R'\(A(:,e)'*W'*W*r));
%x = x + err;

elevation = reshape(x,rows+1,cols+1);
disp(['solution time ' num2str(toc)]);