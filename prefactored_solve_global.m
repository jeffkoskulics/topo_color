function elevation = prefactored_solve_global(u_sn)

%%prefactored_solve computes the elevation based on a prefactored Q-less QR
%%factorization.  Inputs needed are weighting matrix W, triangular matrix
%%R, and the measurement vector u_sn

global W R R_transp A_perm_transp e

tic;
[rows,cols,~] = size(u_sn);

normals(:,1) = reshape(u_sn(:,:,1),rows*cols,1);
normals(:,2) = reshape(u_sn(:,:,2),rows*cols,1);
normals(:,3) = reshape(u_sn(:,:,3),rows*cols,1);

%Identify bad slope measurements and set to zero slope
bad_pixels = lt(normals(:,3),0.5);
normals(bad_pixels,1:2) = 0;  %Assume the slope of bad pixels equals zero
normals(bad_pixels,3)   = 1;  %Zero slope has normal pointing straight up

%Organize measurement vector components
b_x = normals(:,1)./normals(:,3);
b_y = normals(:,2)./normals(:,3);
curv_x = zeros((rows+1)*(cols-1),1); 
curv_y = zeros((rows-1)*(cols+1),1); 
z_mean = 0.63*ones((rows+1)*(cols+1),1);

%Correct nonzero mean
mean_slope_x = mean(b_x);
mean_slope_y = mean(b_y);
b_x = b_x - mean_slope_x; %Assume slope of zero
b_y = b_y - mean_slope_y; 

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