function [surface] = SolveSystem_global(data,u_fp,X,Y,PL,spectra_cube, ...
    flat_surf,fp_res,iteration_vector,z_fp)
% wire: vector for elevation readings [x_pos,y_pos,elevation] from wire/rod
% iteration vector: vector that specifies the desired iteration
%     configuration to run the program in. Specify as slope iterations per
%     elevation iteration. So [1,3,1,2] will run 1 slope -> 1 elevation -> 
%     3 slopes -> 1 elevation -> 1 slope -> 1 elevation ...

%Initialize the variables
u_sn_previous = surface_normals(flat_surf); %Assume flat surface
previous_surface = flat_surf;
%bad_pixels = false(size(u_sn_previous,1),size(u_sn_previous,2)); %Assume no bad pixels

tStart = tic;

for j = 1:length(iteration_vector)
    tic
    
    for i = 1:iteration_vector(j)
        % Retrieve Normals
        u_sn = SlopeRetrieval(data,u_sn_previous, ...
            previous_surface,u_fp,X,Y,PL,spectra_cube);
        u_sn_previous = u_sn;
        
        display(['Slope Iteration #',int2str(i),' complete.'])
        
    end
    
    %Update bad_pixels on the first run
    %bad_pixels = check_normals(u_sn,bad_pixels);
    
    %Retrieve Surface
    elevation = prefactored_solve_global(u_sn, fp_res);

    %Create the surface from the new elevation values
    surface(:,:,1:2) = previous_surface(:,:,1:2); %x and y don't change
    surface(:,:,3) = elevation;

    %Update new surface
    previous_surface = surface;
    u_sn = surface_normals(surface);

    tPassed = num2str(toc);
    display(['Elevation Iteration #',int2str(j),' complete in ',tPassed,' seconds'])
    
end
    
display('Data Analysis Complete')
tPassed = num2str(toc(tStart));
display([10,'Total time for ',int2str(j),' iterations is ',tPassed])