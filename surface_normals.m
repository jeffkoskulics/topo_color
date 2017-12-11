function [u_sn, r_sn] = surface_normals(r_surf)
%%Surface Vectors outputs normalized (u_sn) and unormalized (r_sn) surface
%%vectors
%This calculation returns a matrix of surface normal vectors with one fewer row
%and column than the parent matrix.  Also, this calculation assumes that
%surface normals are well-determined, using only three triangular points.
%u_s(urface)n(ormal) = cross product of vectors going from each point to 
%two of its neighbors (down 1 row/over 1 col).
rows = size(r_surf,1)-1;
cols = size(r_surf,2)-1;

%r_sn = cross([ dx 0 0 ], [0 dy 0]) 
r_sn = cross((r_surf(1:rows,2:cols+1,:)-r_surf(1:rows,1:cols,:)),...
    (r_surf(2:rows+1,1:cols,:)-r_surf(1:rows,1:cols,:)),3);
mag = sqrt(dot(r_sn,r_sn,3));
u_sn(:,:,1) = r_sn(:,:,1) ./ mag; %Normalize
u_sn(:,:,2) = r_sn(:,:,2) ./ mag;
u_sn(:,:,3) = r_sn(:,:,3) ./ mag;

clear mag rows cols

