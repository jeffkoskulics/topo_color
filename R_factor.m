function R_factor(rows,cols,spacing,parent_directory)
% This function produces the Q-less QR factorization of the connectivity
% matrix R and a weighting matrix W

global A_perm_transp e W R R_transp

%Set up the connectivity matrix
disp('Creating matrices...');
[a,b,c,d] = connection(rows,cols);
A = [spacing^-1 * a; spacing^-1 * b; spacing^-2 * c; spacing^-2 * d;...
    speye((rows+1)*(cols+1))];

[num_x_curv ~] = size(c);
[num_y_curv ~] = size(d);

%Standard deviation vector
slope_sigma_x = tan(2.0970*pi/180); %2.0970 degree x-slope standard deviation
slope_sigma_y = tan(2.3477*pi/180); %2.3477 degree y-slope standard deviation
b_x_sigma = slope_sigma_x*ones(rows*cols,1);
b_y_sigma = slope_sigma_y*ones(rows*cols,1);

curv_sigma = 40;
curv_x_sigma = curv_sigma*ones(num_x_curv,1);
curv_y_sigma = curv_sigma*ones(num_y_curv,1);

elev_sigma = 0.05;
z_sigma = elev_sigma*ones((rows+1)*(cols+1),1);

%Full standard deviation vector
w = [b_x_sigma; b_y_sigma; curv_x_sigma; curv_y_sigma; z_sigma].^(-1);

%Weighting matrix C
n = numel(w);
W = sparse(1:n,1:n,w,n,n); %Inverse of diagonal covariance matrix

%Weighted least squares solution
disp('Column reordering...');
A = W*A;
tic;
e = colamd(A);
A_perm_transp = A(:,e)';
clear A
disp(['colamd time ' num2str(toc)]);
disp('Factoring...');
tic 
R = qr(A_perm_transp',0);
disp(['QR factorization time ' num2str(toc)]);
disp('Creating R transpose')
R_transp = R';

%Save the R W e and C factors in a .mat file (several GB)
disp('Saving factor...'); tic;
save([parent_directory, 'R_factors.mat'], 'R', 'R_transp','W','e','A_perm_transp','-v7.3');
disp(['Saving complete. Saving time:  ' num2str(toc)]);

clear global A_perm_transp e W R R_transp