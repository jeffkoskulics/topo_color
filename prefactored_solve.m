function elevation = prefactored_solve_global(u_sn, fp_res)

%%prefactored_solve computes the elevation based on a prefactored Q-less QR
%%factorization.  Inputs needed are weighting matrix W, triangular matrix
%%R, and the measurement vector u_sn

global W R A e

tic;
[rows,cols,~] = size(u_sn);
normals(:,1) = reshape(u_sn(:,:,1),rows*cols,1);
normals(:,2) = reshape(u_sn(:,:,2),rows*cols,1);
normals(:,3) = reshape(u_sn(:,:,3),rows*cols,1);
bad_pixels = lt(normals(:,3),0.5);
normals(bad_pixels,1:2) = 0;  %Assume the slope of bad pixels equals zero
normals(bad_pixels,3)   = 1;  %Zero slope has normal pointing straight up

b_x = normals(:,1)./normals(:,3);
b_y = normals(:,2)./normals(:,3);
curv_x = zeros((rows+1)*(cols-1),1); 
curv_y = zeros((rows-1)*(cols+1),1); 
z_mean = 0.63*ones((rows+1)*(cols+1),1);

b = [b_x; b_y; curv_x; curv_y; z_mean];

tic;
%x = (W*A)\(W*b);
x(e,1) = R\(R'\(A(:,e)'*W'*W*b));
%r = W*b - W*A*x;
%err(e,1) = R\(R'\(A(:,e)'*W'*W*r));

%x = x + err;

elevation = reshape(x,rows+1,cols+1);
disp(['solution time ' num2str(toc)]);