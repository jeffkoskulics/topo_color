function [row_slope col_slope row_curv col_curv jacobian] = connection(rows,cols)

%Universal connection matrix (to test memory requirements for
%factorization)
%connection_rand replaces the ones and minus ones with random numbers
%connection sets up the sparse matrix, A, that multiplies a state vector
%of elevations, x, to give a measurement vector of slope components and 
%mean elevations, b.  A is sparse and has a block structure that creates 
%a concatenated b vector 
%y-slopes are concatenated vertically with x-slopes and mean elevation
%points:
%
% A = [Ay; Ax; I] and b = [by; bx; mean_elevation]
%
%Here, the mean elevation is a vector whose number of elements matches
%the number of elevation vertices in the state vector
%
%rows = the number of rows of a data image
%cols = the number of columns in a data image

%Setup a matrix of measurement indices (subscripts)
row_index = repmat((1:rows)',1,cols);
col_index = repmat(1:cols,rows,1);
measurement_index = sub2ind([rows cols],row_index,col_index);
measurement_vector = reshape(measurement_index,rows*cols,1);
clear row_index col_index measurement_index

%Setup a matrix of vertex (node) indices (subscripts).  This matrix has a
%size of (rows+1) x (cols+1) 
row_index = repmat((1:rows+1)',1,cols+1);
col_index = repmat(1:cols+1,rows+1,1);
state_index = sub2ind([rows+1 cols+1],row_index,col_index);
clear row_index col_index

%Setup a vector containing the vertex indices for the upper left, lower
%left, and upper right vertices
upper_left = reshape(state_index(1:rows,1:cols),rows*cols,1);
lower_left = reshape(state_index(2:rows+1,1:cols),rows*cols,1);
upper_right = reshape(state_index(1:rows,2:cols+1),rows*cols,1);
num_vertices = numel(state_index);

%Sparse matrix for the y-slope components 
A_up_left = sparse(measurement_vector,upper_left,ones(size(measurement_vector)),...
        numel(measurement_vector),num_vertices);
A_low_left = sparse(measurement_vector,lower_left,-ones(size(measurement_vector)),...
        numel(measurement_vector),num_vertices);
A_up_right = sparse(measurement_vector,upper_right,-ones(size(measurement_vector)),...
        numel(measurement_vector),num_vertices);
col_slope = A_up_left + A_low_left;
row_slope = A_up_left + A_up_right;
jacobian = A_up_left + A_low_left + A_up_right;
%A = [A;B];
clear A_up_left A_low_left A_up_right

%Sparse matrix for elevations (identity matrix)
%A_elev = speye(num_vertices,num_vertices);

%Combine y- and x- matrices
%A = [A; A_elev];

%2nd order finite difference equations for the horizontal curvature

%Setup a vector containing the vertices for 2nd order fininte difference
%equations for horizontal curvature.  It will look at three stripes spanning the
%vertical dimension
hor_left = reshape(state_index(:,1:cols-1),(rows+1)*(cols-1),1);
hor_mid  = reshape(state_index(:,2:cols)  ,(rows+1)*(cols-1),1);
hor_rite = reshape(state_index(:,3:cols+1),(rows+1)*(cols-1),1);
meas = 1:numel(hor_left);
Hor_Curv_1 = sparse(meas,hor_left,ones(size(hor_left)),numel(meas),num_vertices);
Hor_Curv_2 = sparse(meas,hor_mid,-2*ones(size(hor_left)),numel(meas),num_vertices);
Hor_Curv_3 = sparse(meas,hor_rite,ones(size(hor_left)),numel(meas),num_vertices);
row_curv = Hor_Curv_1 + Hor_Curv_2 + Hor_Curv_3;
clear Hor_Curv_1 Hor_Curv_2 Hor_Curv_3 hor_left hor_mid hor_rite

%Setup a vector containing the vertices for 2nd order finite differences
%for vertical curvature.  It will look at three bands spanning the
%horizontal dimension

vert_top = reshape(state_index(1:rows-1,:),(rows-1)*(cols+1),1);
vert_mid = reshape(state_index(2:rows  ,:),(rows-1)*(cols+1),1);
vert_bot = reshape(state_index(3:rows+1,:),(rows-1)*(cols+1),1);
meas = 1:numel(vert_top);
Vert_Curv_1 = sparse(meas,vert_top,ones(size(vert_top)),numel(meas),num_vertices);
Vert_Curv_2 = sparse(meas,vert_mid,-2*ones(size(vert_top)),numel(meas),num_vertices);
Vert_Curv_3 = sparse(meas,vert_bot,ones(size(vert_top)),numel(meas),num_vertices);
col_curv = Vert_Curv_1 + Vert_Curv_2 + Vert_Curv_3;
clear Vert_Curv_1 Vert_Curv_2 Vert_Curv_3


