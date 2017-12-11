function bad_pixels = check_normals(u_sn,bad_pixels)
% check_normals checks for surface normals that appear to be
% invalid. Validity is checked by looking for locations where
% the z component of the surface normal is less than .5

bad_pixels = bad_pixels | lt(u_sn(:,:,3),.5);