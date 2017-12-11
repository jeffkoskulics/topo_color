%ConvergenceTests tests the speed of convergence of our model under
%different conditions

%{ 
HOW TO CONFIGURE THIS SCRIPT
The main purpose of this script is to test the speed of convergence
under a number of different iteration conditions. Since the surface is
solved using a two part model, it may be advantageous to run multiple
iterations of the first step before ever running the second.

Since you can only run at most 1 elevation retrieval, regardless of the
number of slope retrievals, your input into the script is in the units:
    slope retrievals per elevation retrieval.

Since at least one slope retrieval should be run for every elevation
retrieval you must NEVER specify a 0!

You should create a variable named "run" as a STRING of
numbers. Each new row should be a new string using the format:

run{1} = '#########';
run{2} = '#########';
...

A sample run is as follows:   run{1} = '555113';
    Which will run 
        5 Slopes > 1 Elevation > 5 Slopes > 1 Elevation >
        5 Slopes > 1 Elevation > 1 Slope  > 1 Elevation >
        1 Slope  > 1 Elevation > 3 Slopes > 1 Elevation

All outputs will be saved in a  subdirectory of the directory of the 
script. Each new row is a run of the test script and will create a new 
directory of outputs.

%}
% Get variable name
data_directory = input([10,'What is the name of the data directory?',10], 's');
output_name = input([10,'Specify a unique name for the data run.',10], 's');

% Load in the first image from directory, along with wire data

% Check to see if directory ends with a '\'
if data_directory(length(data_directory)) ~= '\'
    data_directory = [data_directory,'\'];
end
    
%Load the model constants
load([data_directory,'model_constants.mat'])

%Import the list of files from that directory
display('Importing File List')
File_List = dir([data_directory,'*.tif']);

%Import the instrument data
inst_name = dir([data_directory,'*.dat']);
inst_data = dlmread([data_directory,inst_name(1).name],'\t',1,0);

%Process the inst data to find camera pulses
% The readings are converted from volts to mm
m = 1;
q = 1;
while q <= length(inst_data)-30
    if inst_data(q,18) > 5
        elevation_A(m,1) = (mean(inst_data(q-5:q+15,14))*3.243)+633;
        elevation_B(m,1) = (mean(inst_data(q-5:q+15,15))*2.294)+633;
        elevation_C(m,1) = (mean(inst_data(q-5:q+15,16))*2.985)+633;
        elevation_D(m,1) = (mean(inst_data(q-5:q+15,17))*3.327)+633;
        m = m + 1;
        q = q + 30; %There are about 68 points between pulses
                    % and 7 points per pulse.
    end
    
    q = q + 1;
end

%Load in the second image and the second line of elevation data
wire = zeros(4,3);
wire = [45 53 elevation_A(2); 989 57 elevation_B(2); 40 1230 elevation_D(2); 1002 1239 elevation_C(2)];
wire(:,3) = wire(:,3)*0.001;

surfsolvedata = load_data(data_directory,File_List(2).name);
surfsolvedata = double(demosaic(uint16((double(surfsolvedata) - offset)),'gbrg')).*gain;


if ~exist([data_directory,'convergence_tests'],'dir')
        mkdir([data_directory,'convergence_tests']);
end

save_dir = [[data_directory,'convergence_tests\'],output_name];

if ~exist(save_dir,'dir')
        mkdir(save_dir);
end

%Pull system constants to workspace
system_constants

%Check to see if data variable already demosaiced
% if so, assume the flat_cap was already subtracted
imagesc(uint16(surfsolvedata));
title('Raw Data Image')
saveas(gcf,[save_dir,'\0 - DataImage.tif'],'tiff');
close(gcf);

results = cell(size(run,2),2);

%names of each column
col_names{1,1} = 'Iteration Number';
col_names{1,2} = 'Total Time (sec)';
col_names{1,3} = 'SN X Average';
col_names{1,4} = 'SN Y Average';
col_names{1,5} = 'SN X RMS';
col_names{1,6} = 'SN Y RMS';
col_names{1,7} = 'SN X Std Dev';
col_names{1,8} = 'SN Y Std Dev';
col_names{1,9} = 'Data Red RMS';
col_names{1,10} = 'Data Green RMS';
col_names{1,11} = 'Model Red Average';
col_names{1,12} = 'Model Green Average';
col_names{1,13} = 'Model Red RMS';
col_names{1,14} = 'Model Green RMS';
col_names{1,15} = 'Model Red Std Dev';
col_names{1,16} = 'Model Green Std Dev';
col_names{1,17} = 'Elevation Average';
col_names{1,18} = 'Elevation RMS';
col_names{1,19} = 'Elevation Std Dev';

%Calculate flat model image to data image residual
bad_pixels = false(1024,1280);
flat_normals = surface_normals(flat_surf);
flat_model = forward_model(flat_normals,flat_surf,u_fp,1,X,Y,PL,spectra_cube,bad_pixels);
flat_model_data_rms = check_residual(flat_model,cat(3,surfsolvedata(:,:,1).*~bad_pixels,surfsolvedata(:,:,2).*~bad_pixels,surfsolvedata(:,:,3)));
clear flat_normals flat_model

%Each row in "run" is a new run of this script
for i = 1:size(run,2)
    
    %Number of rows used when recording data
    row_count = 1;
    slope_iter = 0;
    elev_iter = 0;
    
    %Input the name for this run
    results{i,1} = run{i};
    
    %Pre allocate for speed
    size_of_run = 0;
    for p = 1:length(run{i})
        size_of_run = size_of_run + str2double(run{i}(p)) + 1;
    end
    results{i,2} = zeros(size_of_run,19);
    
    %{
    %Output directory for this run
    subdir = [run{i},'\'];
    if ~exist([save_dir,'\',subdir],'dir')
        mkdir([save_dir,'\',subdir]);
    end
    %}
    
    %Assume a flat surface to start
    u_sn_previous = surface_normals(flat_surf);
    previous_surface = flat_surf;
    previous_test_image = surfsolvedata;
    bad_pixels = false(size(u_sn_previous,1),size(u_sn_previous,2)); %Assume no bad pixels
    
    display([10,'***Starting Test Case: ',run{i},'***'])
    
    tic
    %The length of the string in "run" gives the number of elevation
    for j = 1:size(run{i},2)
        num_slopes = str2double(run{i}(j));
        
        %The number of slopes per iteration is given by the single integer
        % value in the string
        for k = 1:num_slopes
            
            u_sn = SlopeRetrieval(surfsolvedata,u_sn_previous,previous_surface,u_fp,X,Y,PL,spectra_cube,bad_pixels);
            
            tPassed = toc;
            slope_iter = slope_iter + 1;
            
            results{i,2}(row_count,1) = 0; %denotes slope iteration
            results{i,2}(row_count,2) = tPassed; %total time passed since start of this run
            
            %Calculate residual x and y for surface normals
            [rms,std_dev,average] = check_residual(u_sn,u_sn_previous);
            results{i,2}(row_count,3) = average(1); %sn_x_average
            results{i,2}(row_count,4) = average(2); %sn_y_average
            results{i,2}(row_count,5) = rms(1); %sn_x_rms
            results{i,2}(row_count,6) = rms(2); %sn_y_rms
            results{i,2}(row_count,7) = std_dev(1); %sn_x_std_dev
            results{i,2}(row_count,8) = std_dev(2); %sn_y_std_dev
            
            bad_pixels = check_normals(u_sn,bad_pixels);
           
            %Create an image and calc residual for that
            test_image = forward_model(u_sn,previous_surface,u_fp,0,X,Y,PL,spectra_cube,bad_pixels);
            figure
            imagesc(uint16(test_image));
            title(['Run',run{i},' S',int2str(slope_iter),' E',int2str(elev_iter)])
            saveas(gcf,[save_dir,'\1 - Run',run{i},'_',int2str(row_count),'_S',int2str(slope_iter),'_E',int2str(elev_iter),'.tif'],'tiff');
            close(gcf);
            figure
            imagesc(uint16(2^6*(cat(3,surfsolvedata(:,:,1).*~bad_pixels,surfsolvedata(:,:,2).*~bad_pixels,surfsolvedata(:,:,3)) - test_image)));
            title(['Run',run{i},' S',int2str(slope_iter),' E',int2str(elev_iter),' Differences'])
            saveas(gcf,[save_dir,'\2 - Run',run{i},'_',int2str(row_count),'_S',int2str(slope_iter),'_E',int2str(elev_iter),'_Differences','.tif'],'tiff');
            close(gcf);
            
            %Between real data image and newest guess
            [rms] = check_residual(test_image,cat(3,surfsolvedata(:,:,1).*~bad_pixels,surfsolvedata(:,:,2).*~bad_pixels,surfsolvedata(:,:,3)));
            results{i,2}(row_count,9) = rms(1); %data_red_residual
            results{i,2}(row_count,10) = rms(2); %data_green_residual
            
            %Between newest guess and last data image
            [rms,std_dev,average] = check_residual(test_image,previous_test_image);
            results{i,2}(row_count,11) = average(1); %model_red_average
            results{i,2}(row_count,12) = average(2); %model_green_average
            results{i,2}(row_count,13) = rms(1); %model_red_rms
            results{i,2}(row_count,14) = rms(2); %model_green_rms
            results{i,2}(row_count,15) = std_dev(1); %model_red_std_dev
            results{i,2}(row_count,16) = std_dev(2); %model_green_std_dev
            
            if j == 1
                %first run, no elevation results
                results{i,2}(row_count,17) = 0; %elevation_average
                results{i,2}(row_count,18) = 0; %elevation_rms
                results{i,2}(row_count,19) = 0; %elevation_std_dev
            else
                %Set equal to previous elevation's results
                results{i,2}(row_count,17) = results{i,2}(row_count-1,17); %elevation_average
                results{i,2}(row_count,18) = results{i,2}(row_count-1,18); %elevation_rms
                results{i,2}(row_count,19) = results{i,2}(row_count-1,19); %elevation_std_dev
            end
            
            previous_test_image = test_image;
            u_sn_previous = u_sn;
            row_count = row_count + 1;
            
            %{
            figure
            nnz(bad_pixels)
            imagesc(bad_pixels)
            %}
            
            display(['Slope Iteration: ',int2str(slope_iter),' Completed.'])
        end
        
        slope_iter = 0;
        
        [surface,u_sn,bad_pixels] = ElevationRetrieval(u_sn,previous_surface,fp_res,wire,bad_pixels,z_fp);
        
        elev_iter = elev_iter + 1;
        tPassed = toc;
        results{i,2}(row_count,1) = 1; %elevation iteration
        results{i,2}(row_count,2) = tPassed; %overall time passed
        
        %Calculate residual elevation
        [rms,std_dev,average] = check_residual(surface(:,:,3),previous_surface(:,:,3));
        results{i,2}(row_count,17) = average(1); %elevation_average
        results{i,2}(row_count,18) = rms(1); %elevation_rms
        results{i,2}(row_count,19) = std_dev(1); %elevation_std_dev
        
        %Create an image and calc residual for that
        figure
        test_image = forward_model(u_sn,previous_surface,u_fp,0,X,Y,PL,spectra_cube,bad_pixels);
        imagesc(uint16(test_image));
        title(['Run',run{i},' S',int2str(slope_iter),' E',int2str(elev_iter)])
        saveas(gcf,[save_dir,'\3 - Run',run{i},'_',int2str(row_count),'_E',int2str(elev_iter),'.tif'],'tiff');
        close(gcf);
        figure
        imagesc(2^6*(uint16(cat(3,surfsolvedata(:,:,1).*~bad_pixels,surfsolvedata(:,:,2).*~bad_pixels,surfsolvedata(:,:,3)) - test_image)));
        title(['Run',run{i},' S',int2str(slope_iter),' E',int2str(elev_iter),' Differences'])
        saveas(gcf,[save_dir,'\4 - Run',run{i},'_',int2str(row_count),'_E',int2str(elev_iter),'_Differences','.tif'],'tiff');
        close(gcf);
        
        %Between real data image and newest guess
        [rms] = check_residual(test_image,cat(3,surfsolvedata(:,:,1).*~bad_pixels,surfsolvedata(:,:,2).*~bad_pixels,surfsolvedata(:,:,3)));
        results{i,2}(row_count,9) = rms(1); %data_red_residual
        results{i,2}(row_count,10) = rms(2); %data_green_residual
            
        %Between newest guess and last data image
        [rms,std_dev,average] = check_residual(test_image,previous_test_image);
        results{i,2}(row_count,11) = average(1); %model_red_average
        results{i,2}(row_count,12) = average(2); %model_green_average
        results{i,2}(row_count,13) = rms(1); %model_red_rms
        results{i,2}(row_count,14) = rms(2); %model_green_rms
        results{i,2}(row_count,15) = std_dev(1); %model_red_std_dev
        results{i,2}(row_count,16) = std_dev(2); %model_green_std_dev
        
        %Update surface normals
        [rms,std_dev,average] = check_residual(u_sn,u_sn_previous);
        results{i,2}(row_count,3) = average(1); %sn_x_average
        results{i,2}(row_count,4) = average(2); %sn_y_average
        results{i,2}(row_count,5) = rms(1); %sn_x_rms
        results{i,2}(row_count,6) = rms(2); %sn_y_rms
        results{i,2}(row_count,7) = std_dev(1); %sn_x_std_dev
        results{i,2}(row_count,8) = std_dev(2); %sn_y_std_dev
        
        
        %Plot current surface
        plot_surface(surface);
        title(['Run',run{i},' E',int2str(elev_iter),' Surface'])
        saveas(gcf,[save_dir,'\5 - Run',run{i},' E',int2str(elev_iter),' Surface','.tif'],'tiff');
        close(gcf);
        
        %Plot differences between current and previous surface
        if elev_iter ~= 1
            plot_surface_diff(surface,previous_surface);
            title(['Run',run{i},' E',int2str(elev_iter),' - E',int2str(elev_iter-1),' Difference'])
            saveas(gcf,[save_dir,'\6 - Run',run{i},' E',int2str(elev_iter),' - E',int2str(elev_iter-1),' Difference','.tif'],'tiff');
            close(gcf);
        end
        
        previous_test_image = test_image;
        previous_surface = surface;
        u_sn_previous = u_sn;
        row_count = row_count + 1;
            
        
        %{
        figure
        nnz(bad_pixels)
        imagesc(bad_pixels)
        %}
        
        display(['Elevation Iteration: ',int2str(elev_iter),' Completed.'])
    end
    
    %%%  Make and save plots of aggregate data  %%
    
    %Parse run name
    temp_string = '';
    for q = 1:size(run{i},2)-1
        temp_string = [temp_string,run{i}(q),','];
    end
    temp_string = [temp_string,run{i}(size(run{i},2))];
    
    %Plot surface normal RMS
    figure;
    plot(results{i,2}(:,2),results{i,2}(:,6:7))
    legend('Surface Normal X','Surface Normal Y')
    title(['Run: ',temp_string,' | RMS Difference of SN Iterations'])
    hold on
    xlabel('Time (seconds)')
    for q = 1:size(results{i,2},1)
        if results{i,2}(q,1) == 1
            plot(results{i,2}(q,2),results{i,2}(q,6),'Marker','x','MarkerFaceColor','k','MarkerSize',10)
            plot(results{i,2}(q,2),results{i,2}(q,7),'Marker','x','MarkerFaceColor','k','MarkerSize',10)
        end
    end
    hold off
    saveas(gcf,[save_dir,'\Run',run{i},'_SN_RMS.fig']);
    close(gcf);
    
    %Plot Data Image to Model Image RMS Difference
    figure;
    plot([0;results{i,2}(:,2)],[flat_model_data_rms(1:2);results{i,2}(:,9:10)]) %adding in the zero point
    legend('Red Channel','Green Channel')
    title(['Run: ',temp_string,' | RMS Difference of Data Image from Model Image'])
    hold on
    xlabel('Time (seconds)')
    for q = 1:size(results{i,2},1)
        if results{i,2}(q,1) == 1
            plot(results{i,2}(q,2),results{i,2}(q,9),'Marker','x','MarkerFaceColor','k','MarkerSize',10)
            plot(results{i,2}(q,2),results{i,2}(q,10),'Marker','x','MarkerFaceColor','k','MarkerSize',10)
        end
    end
    hold off
    saveas(gcf,[save_dir,'\Run',run{i},'_DataAndModel_RMS.fig']);
    close(gcf);
    
    %Plot Model Image to Model Image RMS Difference
    figure;
    plot(results{i,2}(:,2),results{i,2}(:,13:14))
    legend('Red Channel','Green Channel')
    title(['Run: ',temp_string,' | RMS Difference of Model Image Iterations'])
    hold on
    xlabel('Time (seconds)')
    for q = 1:size(results{i,2},1)
        if results{i,2}(q,1) == 1
            plot(results{i,2}(q,2),results{i,2}(q,13),'Marker','x','MarkerFaceColor','k','MarkerSize',10)
            plot(results{i,2}(q,2),results{i,2}(q,14),'Marker','x','MarkerFaceColor','k','MarkerSize',10)
        end
    end
    hold off
    saveas(gcf,[save_dir,'\Run',run{i},'_ModelAndModel_RMS.fig']);
    close(gcf);
    
    %Plot Elevation RMS Difference
    figure;
    elev_start = 1;
    for q = 1:size(results{i,2},1)
        if results{i,2}(q,1) == 1
            break
        end
        elev_start = elev_start + 1;
    end
    %leave off first points where there is no elevation
    plot(results{i,2}(elev_start:size(results{i,2},1),2),results{i,2}(elev_start:size(results{i,2},1),18))
    legend('Elevation (meters)')
    title(['Run: ',temp_string,' | RMS Difference of Elevation Iterations'])
    hold on
    xlabel('Time (seconds)')
    for q = 1:size(results{i,2},1)
        if results{i,2}(q,1) == 1
            plot(results{i,2}(q,2),results{i,2}(q,18),'Marker','x','MarkerFaceColor','k','MarkerSize',10)
        end
    end
    hold off
    saveas(gcf,[save_dir,'\Run',run{i},'_Elevation_RMS.fig']);
    close(gcf);
    
    %Save final elevation
    elevation = surface(:,:,3);
    save([save_dir,'\Run',run{i},'_elevation.mat'],'elevation')
end

%Save the analysis matrix
save([save_dir,'\Analysis_Results.mat'],'results','col_names')

clear X Y cam_pos flat_fitting flat_image flat_surf fp_res ...
        lens_cap nadir r_ap u_fp z_aperture wire clear_const z_fp ...
        average bad_pixels elev_iter i j k num_slopes rms slope_iter ...
        std_dev var_name surfsolvedata test_image u_sn_previous output_name ...
        p previous_surface previous_test_image row_count save_dir size_of_run ...
        subdir tPassed





