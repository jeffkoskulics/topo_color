function elevation = solve_limited_data4(u_sn,previous_elevation,...
    bad_pixels,wire,fp_res)
%% solve_limited_data finds a topography that best fits surface normal unit 
%  vector data, excluding any (bad) pixels that are flagged as unreliable
%  slopes.  Unconnected vertices are set equal to the elevation from the
%  previous iteration
%  
%  this 4th version gives each pixel in the elevation box its own row, as
%  opposed to putting all pixels in a row to average together
%
%  A - sparse connectivity matrix
%  u_sn - measurement matrix (including unreliable surface normals)
%  previous_elevation - previous iteration's assumed surface
%  wire = [rows cols elevs], capacitive wire rows, columns and elevations
%  in array form for the image beign solved
%  bad_pixels - binary array (1's indicate bad pixels)

%Create the connectivity matrix
[rows cols k] = size(u_sn); 
box = 10;
A = connection4(rows,cols,wire(:,1),wire(:,2),box);
clear k
%Create a vector of wire measurements repeated for every pixel in each box
N = 2*box;
elevs = zeros(4*N^2,1);
for i = 1:size(wire,1)
    elevs(N^2*(i-1)+1:N^2*i) = repmat(wire(i,3),N^2,1);
end
%Create a vector of slopes in three columns ([x y z] components)
normals = reshape(u_sn,rows*cols,3);
%Scale the components into differences of elevation
normals(:,1) = fp_res*normals(:,1)./normals(:,3); %Check to make sure 
normals(:,2) = fp_res*normals(:,2)./normals(:,3); %connection matches slopes
%Patch the slopes vectors [y-component; x-component; elevations]
b = [normals(:,2);normals(:,1);elevs]; 

num_pixels = numel(bad_pixels); 
%Linear index of bad pixels
bad_rows = find(bad_pixels); 
%Both x- and y-normal components must be removed
bad_rows = [bad_rows; bad_rows + num_pixels]; 

%Remove all rows connected to bad data
A(bad_rows,:) = []; 
%Remove bad data from measurement vector
b(bad_rows) = [];   

%Columns with zero Euclidean length (excluding the elevation rows at the
% end) are considered vertices that are not connected to any other vertex.  
good_slopes = size(A,1) - size(elevs,1);
bad_columns = find(lt(dot(A(1:good_slopes,:),A(1:good_slopes,:)),1)); 

%Save the index values of connected elevations
good_columns = 1:size(A,2);
good_columns(bad_columns) = [];     

%Delete all columns of unconnected vertices
A(:,bad_columns) = []; 

%Weight the rows for elevations according to the number of non-zero entries
%to correctly calculate the mean elevation (for each wire elevation
%measurement)
%{
for i = 1:size(wire,1)
    good_vertices = nnz(A(good_slopes+i,:));
    A(good_slopes+i,:) = A(good_slopes+i,:)/good_vertices;
    disp('nnz good pixels around elevation')
    [i good_vertices sum(A(good_slopes+i,:)) wire(i,3)]
end
disp('Sum of rows for slope portion')
sum(abs(sum(A(1:good_slopes,:),2)))
%}

%Display number of remaining columns (subtract rank from the number
%of remaining columns for count of unconnected islands)
disp([num2str(size(A,2)),' columns remaining in A']) 

%Weighting matrix
%Start by constructing a diagonal matrix with elements equal to the
%reciprocal of the estimated measurement standard deviations
W = speye(length(b));

%Slope measurements are assumed to have a standard deviation of one degree
W(1:good_slopes,:) = (fp_res*sin(pi/180)*cos(pi/180))^-1 * W(1:good_slopes,:);
%Elevation measurements are assumed to have a standard deviation of two
%millimeters
W(good_slopes+(1:size(elevs)),:) = 0.002^-1*W(good_slopes+[1:size(elevs,1)]',:);

%Set the elements equal to variance (the standard deviation squared)
W = W.^2;

%Multiply the weighting matrix by both sides of normal equation:
% W*A*x = W*b;
A = W*A;
b = W*b;
%b(length(b)-10:length(b))

%Least squares solution
z = A\b;    

elevation = zeros(numel(previous_elevation),1);

%Reconstruct the solution vector with zeros for all unknown elements;
elevation(good_columns) = z;  

%******previous_elevation(bad_columns); 
%Set unknown elements equal to previous iteration's elevation
elevation(bad_columns) = 0;  

%Add the unconnected corner pixel (assume it's equal to the elevation in
%the vertex located in the row directly above, in the same column)
elevation(numel(elevation)) = elevation(numel(elevation)-1); 
elevation = reshape(elevation,1025,1281);
%figure; plot(elevation(45,:))
%{
disp('nnz bad_pixels')
nnz(bad_pixels)
disp('max min')
max(max(elevation))
min(min(elevation))
figure; plot(elevation(45,:))
%}