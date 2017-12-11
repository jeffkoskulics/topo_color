function [spectra_cube,X,Y,PL] = SpectrumCalc3(directory,z_fp,fp_res,nadir,r_ap)
%% SpectrumCalc takes in the system properies and creates model to be used in the forward model
% This script takes the input directory and pulls the spectrum data, path
% absorbance data, and camera sensativity data to create a proper model
% cube of system response with respect to position on the backlight and
% path length through the water. The system then scales this cube to match
% the response we see in a flat surface image by running the forward model
% retrieval on unscaled data and correcting the cube before output.

%% %%%%% EDIT SYSTEM CONSTANTS HERE %%%%%
% COORDINATE SYSTEM SUMMARY
% The system is currently set up using a nadir coordinate system in meters. 
% The nadir point is assumed to be the position of (x,y) = (0,0). Positive 
% x runs downtank from the left-hand to right-hand side of the image.
% Positive y runs crosstank from the top to the bottom of the image.
% Basically increasing column number = increaing x, and increasing row
% number = increasing y.

% DIRECTORY and DATA STORAGE
% The spectrum data should be organized into single directory specified by 
% the 'directory' input variable. The 'directory' should have subdirectories
% labeled R#C# corresponding to the row and column of data collection.
% These subdirs should then contain several .dat files with copies of
% several data runs from that same location. Run 'AverageSpectra.m' on the
% parent directory to average these .dat files and convert to .mat files.
% This script expects the .mat outputs from 'AverageSpectra.m'. The parent
% directory should also contain a matrix called 'AbsorptionSensativity.mat'
% which contains the absorbtion and camera sensativity curves.

% SPECTRA LOCATIONS
% Input these locations as vectors running from left to right for x and from
% top to bottom for y. NOTE: the current system has a bad row and column 
% around the perimeter of the image. These are ignored in the loops to 
% load and clear the Spectra data. You may need to change this.
X_locations = [-.2886 -.1366 -.04445 .0476 .1397 .2318];
Y_locations = [-.3023 -.216 -.1295 -.0432 .0432 .1295 .216 .3023];
rows = 8; %Number of row measurements 
cols = 10; %Number of column measurements

% FILTER BOUNDARIES
% These are the measurements of the filter boundaries out to the black bars
% in meters.
filter_top = -.3667;
filter_bottom = .3699;
filter_left = -.4318;
filter_right = .4318;

% WAVELENGTH RANGE
% These should be given in terms of the row numbers in the data files. 
% Since the camera is only sensative in a certain range, and the spectrum 
% range is larger than that, we only need to pull in the visible wavelengths.
start_row = 258;
end_row = 1743;

% Pixel resolution of the flat image
flat_rows = 1024;
flat_cols = 1280;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% Prepare the inputs

%directory should be parent directory of spectrum subdirectories
if directory(length(directory)) ~= '\'
    directory = [directory,'\'];
end

load([directory,'AbsorptionSensitivity.mat'])

%% Fit the Spectra data by x,y, and wavelength

%Load the good spectra data. (Outside row and column are bad)
Spectra = cell(rows,cols);
for i = 2:rows-1
    for j = 2:cols-1
        Spectra{i,j} = dlmread([directory,'R',num2str(i),'C',num2str(j),'.txt'],char(9),[17 0 3664 1]);
    end
end

%X and Y locations of the spectrum readings in the nadir coordinate system
X_vec = repmat(X_locations,length(Y_locations),1);
Y_vec = repmat(Y_locations',1,size(X_vec,2));
X_vec = reshape(X_vec,size(X_vec,1)*size(X_vec,2),1);
Y_vec = reshape(Y_vec,size(Y_vec,1)*size(Y_vec,2),1);

m = 1;
for z = start_row:end_row %The usable wavelengths

    %Pull data for this wavelength only
    q = 1;
    for j = 2:cols-1
        for i = 2:rows-1
            spec(q,1) = Spectra{i,j}(z,2);
            q = q + 1;
        end
    end
    
    % Fit the values
    spectra_fitting(m) = polyfitn([X_vec Y_vec],spec,3);
    m = m + 1;

end

%% Define x,y locations to evaluate polyval at

%Use filter boundaries in nadir coordinate system.
x_start = filter_left;
x_end = filter_right;
x_numpoints = 100;


y_start = filter_top;
y_end = filter_bottom;
y_numpoints = 100;


x_pts = x_numpoints - 1;
y_pts = y_numpoints - 1;
X = x_start:(x_end-x_start)/x_pts:x_end;
Y = y_start:(y_end-y_start)/y_pts:y_end;

poly_val_locations(:,:,1) = repmat(X,length(Y),1);
poly_val_locations(:,:,2) = repmat(Y',1,length(X));

for i = 1:length(spectra_fitting)
    fitted_spectrum(i,:) = polyvaln(spectra_fitting(i),reshape(poly_val_locations,length(X)*length(Y),2));
end
fitted_spectrum = reshape(fitted_spectrum,size(fitted_spectrum,1),x_numpoints,y_numpoints);

%% Create data cube of Response in terms of elevation, x, and y
PL = .50:.05:1.5; %Range of path lengths to sample over
spectra_cube = zeros(x_numpoints,y_numpoints,length(PL),3);
for z = 1:length(PL);
    for x = 1:x_numpoints
        for y = 1:y_numpoints
            spectra_cube(y,x,z,1) = sum(fitted_spectrum(:,y,x).*10.^-(absorption*PL(z)).*red_sens);
            spectra_cube(y,x,z,2) = sum(fitted_spectrum(:,y,x).*10.^-(absorption*PL(z)).*green_sens);
            spectra_cube(y,x,z,3) = sum(fitted_spectrum(:,y,x).*10.^-(absorption*PL(z)).*blue_sens);
        end
    end
end

%% Create grid for data cube interpolation
X_out = repmat(X,[y_numpoints 1 length(PL)]);
Y_out = repmat(Y',[1 x_numpoints length(PL)]);
PL_out(1,1,:) = PL;
PL_out = repmat(PL_out,[x_numpoints y_numpoints 1]);

%{
Scale to full 16 bit range, per channel
maximum = max(max(max(spectra_cube)));
spectra_cube(:,:,:,1) = spectra_cube(:,:,:,1)*(65536/maximum(1));
spectra_cube(:,:,:,2) = spectra_cube(:,:,:,2)*(65536/maximum(2));
spectra_cube(:,:,:,3) = spectra_cube(:,:,:,3)*(65536/maximum(3));
%}

%% Clear all of the raw spectra data
for i = 2:rows-1
    for j = 2:cols-1
        eval(['clear spec R',num2str(i),'C',num2str(j),'Spectra'])
    end
end

%% Create a flat surface and flat surface normals at focal plane elevation
% Focal plane position for each pixel
%r_fp(i,j,:) = [x y z], is calculated using the following logic
%x = (j-nadir(2))*fp_res , y = (i-nadir(1))*fp_res , z = z_fp
r_fp = zeros(flat_rows,flat_cols,3);
r_fp(:,:,1) = repmat((((1:1:flat_cols)-nadir(2))*fp_res),flat_rows,1);
r_fp(:,:,2) = repmat((((1:1:flat_rows)-nadir(1))*fp_res)',1,flat_cols);
r_fp(:,:,3) = repmat(z_fp,flat_rows,flat_cols);

% Incident Vector 
%located on focal plane pointing along the line from the 
%camera aperature to the pixel location.
%u_fp(i,j,:) = fp_positions(i,j,:) - [0,0,z_aperture];
u_fp = zeros(flat_rows,flat_cols,3);
u_fp(:,:,1) = r_fp(:,:,1) - r_ap(1);
u_fp(:,:,2) = r_fp(:,:,2) - r_ap(2);
u_fp(:,:,3) = r_fp(:,:,3) - r_ap(3);
mag = sqrt(dot(u_fp,u_fp,3));
u_fp(:,:,1) = u_fp(:,:,1) ./ mag; %Normalize
u_fp(:,:,2) = u_fp(:,:,2) ./ mag;
u_fp(:,:,3) = u_fp(:,:,3) ./ mag;

%Create flat surface
flat_surf = zeros(flat_rows+1,flat_cols+1,3);
flat_surf(:,:,1) = repmat((((1:1:flat_cols+1)-nadir(2))*fp_res),flat_rows+1,1);
flat_surf(:,:,2) = repmat((((1:1:flat_rows+1)-nadir(1))*fp_res)',1,flat_cols+1);
flat_surf(:,:,3) = repmat(z_fp,flat_rows+1,flat_cols+1);

u_sn = surface_normals(flat_surf);

%Scale the data cube to the 16 bit range
spectra_cube = spectra_cube * 2^16/max(max(max(max(spectra_cube))));

%create the gain matrix
%gain = CreateGain(flat_image,offset,spectra_cube,X_out,Y_out,PL_out,u_sn,flat_surf,u_fp);

%Clear the old X,Y,PL values for better usability of the output
clear X Y PL
X = X_out;
Y = Y_out;
PL = PL_out;

display([10,'Check the directory you specified for ''model_constants.mat'''])
display('Place this in each data directory, it contains the information necessary to run the analysis')
%save([directory,'model_constants.mat'],'spectra_cube','X','Y','PL','flat_surf','u_fp','u_sn','offset','gain')
save([directory,'constants_for_gain3.mat'],'spectra_cube','X_out','Y_out','PL_out','u_sn','flat_surf','u_fp')
