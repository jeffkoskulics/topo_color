function [surface,u_sn] = ElevationRetrieval(u_sn,...
    previous_surface,fp_res,z_fp)
%% %%%% RETRIEVAL STEP 2: ELEVATIONS %%%%
elevation = prefactored_solve(W, R, A, e, u_sn, fp_res);

%Create the surface from the new elevation values
surface = zeros(size(previous_surface,1),size(previous_surface,2),3);
surface(:,:,1:2) = previous_surface(:,:,1:2); %x and y don't change
surface(:,:,3) = elevation;

%Update surface normals based on new surface
u_sn = surface_normals(surface);