function [A] = connection_test(rows,cols,elev_rows,elev_cols)

%connection sets up the sparse matrix connecting vertex elevations with
%slope vector components.
%rows = number of rows of a data image
%cols = number of columns in a data image
%elev_rows = vector of center pixel rows for each capacitive wire elevation
%elev_column = vector of center pixel columns for each capacitive wire elevation

%overwrite elev_rows and elev_cols to assume that 4 wire data will be used
%with the following coordinates (a placeholder to be filled with data
%after testing):
elev_rows = [45 989 40   1002]';
elev_cols = [53 57  1230 1239]';

%Setup a vector of vertex (node) indices
rows = rows + 1;
cols = cols + 1;
row_index = repmat((1:rows)',1,cols);
col_index = repmat(1:cols,rows,1);
node_index = sub2ind([rows cols],row_index,col_index);

%Define a square region around each wire
num_wires = length(elev_rows);
box = 10; %size of square box (in pixels) surrounding each wire
N = 2*box;
wire_node = zeros(N,N,num_wires);
%create a nx2 list containing node indices with corresponding elevations
wire_indices = [0 0]; %a placeholder to concatenate new index values to
for i = 1:num_wires
    %for each wire, find indices of each vertex inside the box
    for j = 1:N
        for k = 1:N
            row = elev_rows(i)-box+j;
            col = elev_cols(i)-box+k;
            wire_node(j,k,i) = sub2ind([rows cols],row,col);
        end
    end
    %reshape into a vector and add to wire_node list
    wire_indices = [wire_indices;reshape(wire_node(:,:,i),numel(wire_node(:,:,i)),1) ...
        repmat(i,N^2,1)];
end
%delete the placeholder
wire_indices(1,:) = []; 

%Delete last row and column so that the size of the vertex index matrix is
%the same size as the edge matrix
node_index(rows,:) = [];
node_index(:,cols) = [];
node_index = reshape(node_index,(rows-1)*(cols-1),1);
node_max = max(node_index);

clear row_index col_index

%Setup a vector of edge (slope vector component) indices
rows = rows - 1;
cols = cols - 1;
row_index = repmat((1:(rows))',1,cols);
col_index = repmat(1:(cols),rows,1);
edge_index = sub2ind([(rows) (cols)],row_index,col_index);
edge_index = reshape(edge_index,(rows)*(cols),1);

points = size(elev_rows,1); 

%Slope vector elements
sparse_row = [edge_index;edge_index;...
              edge_index+(rows*cols);edge_index+(rows*cols)];

%Add rows for the elevation measurements
sparse_row = [sparse_row;max(sparse_row)+wire_indices(:,2)];

%Vertex elevation vector elements     
offset = rows+1;
sparse_col = [node_index;node_index+1;node_index;node_index+offset;...
                wire_indices(:,1)];

%Sparse matrix entries (elements)            
sparse_elements = [ones(size(node_index));-ones(size(node_index));...
    ones(size(node_index));-ones(size(node_index));...
    ones(size(wire_indices,1),1)];

A = sparse(sparse_row,sparse_col,sparse_elements);
