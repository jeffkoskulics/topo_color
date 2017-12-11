% setup_workspace will set the workspace up for you.
% You must load in the flat surface and the lens cap image to
% start the procedure.

display('Setting up Workspace')

% Define system constants
rows = 1024;
cols = 1280;
z_aperture = 3; %meter
cam_pos = [0 0]; %Camera Position [x y] 
z_fp = 0.633; %focal plane elevation, in meters
fp_res = 0.0002898; %in meters on the focal plane
nadir = [1024/2,1280/2]; %pixel location of nadir
r_ap = [cam_pos z_aperture]; %Aperture position

%The rest of the constants are located in model_constants.mat

display('Setup Complete')