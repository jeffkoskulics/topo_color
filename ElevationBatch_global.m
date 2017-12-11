function ElevationBatch_global(parent_dir)
%% ElevationBatch runs the ElevationRetrieval routine on a directory of data
% The data is assumed to be in tif form, and will be output to a
% subdirectory \Output\ with .mat files of the same names as the input
% images. The mat files contain the final surface normals and surface
% elevation matrices.
%Lock-File implementation inspired by: http://www.mathworks.com/matlabcentral/newsreader/view_thread/287163

%% Load and prepare flat surface
clear global
display('Loading constants from data directory')

%%%%INTERMEDIATE OUTPUTS%%%%
%global g_parent_dir g_sub_dir g_file_name g_elev_num

% Generate a process ID
c = clock;
myID = num2str(etime(c,[2012 0 0 0 0 0]));

% Check to see if directory ends with a '\', and convert windows style to
% unix style
[part,remainder] = strtok(parent_dir,'/\');
if isunix
    parent_dir = ['/',part];
else
    parent_dir = part;
end

while ~isempty(part)
    [part,remainder] = strtok(remainder,'/\');
    parent_dir = [parent_dir,'/',part];
end


%Check for output directory existance
if ~exist([parent_dir,'Retrieval/'],'dir')
    display('Output directory does not exist')
    display('Creating output directory')
    mkdir(parent_dir,'Retrieval/');
end

%Load the model and system constants
display('Importing data cube...')
load([parent_dir,'model_constants.mat'])
display('Importing R factors...')
load([parent_dir,'R_factors.mat'])
display('Generating system constants...')
system_constants;
iteration_vector = [2,1];

g_parent_dir = parent_dir;

%% Import Data

%{
 %Instrument data currently ignored in solve_limited_data_5
 %This will require updating to work with the new batch/folder structure
Import the instrument data
inst_name = dir([data_dir,'*.dat']);
inst_data = dlmread([data_dir,inst_name(1).name],'\t',1,0);

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
%}


%% Run the analysis of each data image
display('Begin data image analysis')
complete = false;
previous_data_path = '';

while ~complete
    
    % Lock the job list
    lockSuccess = 0;
    while ~lockSuccess
        % Wait for lock-file to disappear
        fprintf('Checking job-list availability...');
        while exist([parent_dir 'runfile_lock'],'file')
            file_info = dir([parent_dir 'runfile_lock']);
            if (~isempty(file_info) && (now - file_info.datenum > .0035) ) %5 mins ~ .0035 in datenum format
                fprintf('\nLock file more than 5 minutes old, a crash is likely... Deleting lock file.\n');
                delete([parent_dir,'runfile_lock']);
            end
            pause(0.2);
        end
        fprintf(' Available.\n');
        % Attempt to create a lock
        fprintf('Attempting to secure job-list lock file...');
        fid = fopen([parent_dir 'runfile_lock'],'w');
        fprintf(fid,'%s',myID);
        fclose(fid);
        pause(1);
        % Was I successful?
        if strcmp(textread([parent_dir 'runfile_lock'],'%s'),myID)
            lockSuccess = 1;
            fprintf(' Successful.\n');
        else
            fprintf(' Unsuccessful, trying again.\n');
        end
    end

    % Grab the next job
    fprintf('Pulling job from stack...');
    fid = fopen([parent_dir,'runfile.txt'],'r'); 
    z = 1;
    current_line = fgetl(fid);
    [full_path,flag] = strtok(current_line,char(9));

    while (~complete && ~strcmp(flag,[9,'0']))
        if (ischar(current_line))
            z = z + 1;
            current_line = fgetl(fid);
            [full_path,flag] = strtok(current_line,char(9));
        else
            complete = true;
        end
    end
    fclose(fid);
    fprintf(' Job assigned.\n');

    if ~complete
        %Update the run file to show that we are in progress
        fid = fopen([parent_dir,'runfile.txt'],'r+');
        new_line = [full_path,9,'2'];
        % Can change 'loc' to insert data at any line in file. 
        % When loc=2, data will be inserted at line 3 
        loc = z-1; 
        for i = 1:loc 
            temp_line = fgetl(fid);
        end; 
        location = ftell(fid);
        fseek(fid,location,'bof');
        fprintf(fid,'\n'); 
        fseek(fid,-1,'cof'); 

        fprintf(fid,'%s',new_line); 
        fprintf(fid,'\r'); 
        fclose(fid);
    end

    % Remove the job lock
    delete([parent_dir,'runfile_lock']);
    fprintf('Lock file deleted.\n\n');
   
    if ~complete
        % Define the new line
        [full_path,flag] = strtok(current_line,char(9));

        %Parse the full path down to just the file name
        if isunix
            data_path = '/';
        else
            data_path = '';
        end
        token = '';
        file_name = '';
        remain = full_path;
        while ~isempty(remain)
            if ~isempty(token)
                data_path = [data_path,token,'/'];
                file_name = remain(2:length(remain));
            end
            [token,remain] = strtok(remain,'/');
        end
        sub_dir = strrep(data_path,parent_dir,'');
        
        %Check for output directory existance
        if ~exist([parent_dir,'Retrieval/',sub_dir],'dir')
            display('Output directory does not exist')
            display('Creating output directory')
            mkdir([parent_dir,'Retrieval/',sub_dir]);
        end
        
        %%%%INTERMEDIATE OUTPUTS%%%%
        %g_sub_dir = sub_dir;
        %g_file_name = file_name;
        %g_elev_num = 1;

        %Load the image
        if (~strcmp(data_path,previous_data_path))
            display('New directory... loading new gain.mat and offset.mat')
            load([data_path,'gain.mat']);
            load([data_path,'offset.mat']);
        end
        data = load_data(data_path,file_name);
        data = double(demosaic(uint16((double(data) - offset)),'gbrg')).*gain;

        display([10,'Data Image: ',file_name])
        display(['Location: ',data_path,10])

        %Iterate through retrieval
        [surface] = SolveSystem_global(data,u_fp,X,Y,PL,spectra_cube, ...
            flat_surf,iteration_vector);
        surface(:,:,3) = surface(:,:,3) - mean(reshape(surface(:,:,3),...
            numel(surface(:,:,3)),1));

        elevation = single(surface(:,:,3));

        %Save x-y positions if this is the first run
        if (~strcmp(data_path,previous_data_path))
            display('New directory... saving xy_position.')
            xy_position = surface(:,:,1:2);
            save([parent_dir,'Retrieval/',sub_dir,'xy_position.mat'],'xy_position')
        end

        %Save the mat file
        save([parent_dir,'Retrieval/',sub_dir,strtok(file_name,'.'),'.mat'],'elevation')

        % Lock the job list
        lockSuccess = 0;
        while ~lockSuccess
            % Wait for lock-file to disappear
            fprintf('Checking job-list availability...');
            while exist([parent_dir 'runfile_lock'],'file')
                file_info = dir([parent_dir 'runfile_lock']);
                if (~isempty(file_info) && (now - file_info.datenum > .0035) ) %5 mins ~ .0035 in datenum format
                    fprintf('\nLock file more than 5 minutes old, a crash is likely... Deleting lock file.\n');
                    delete([parent_dir,'runfile_lock']);
                end
                pause(0.2);
            end
            fprintf(' Available.\n');
            % Attempt to create a lock
            fprintf('Attempting to secure job-list lock file...');
            fid = fopen([parent_dir 'runfile_lock'],'w');
            fprintf(fid,'%s',myID);
            fclose(fid);
            pause(1);
            % Was I successful?
            if strcmp(textread([parent_dir 'runfile_lock'],'%s'),myID)
                lockSuccess = 1;
                fprintf(' Successful.\n');
            else
                fprintf(' Unsuccessful, trying again.\n');
            end
        end

        %Update the run file to show that we are completed
        fid = fopen([parent_dir,'runfile.txt'],'r+'); 
        new_line = [full_path,9,'1'];

        % Can change 'loc' to insert data at any line in file. 
        % When loc=2, data will be inserted at line 3 
        loc = z-1; 
        for i = 1:loc 
            temp_line = fgetl(fid);
        end; 
        location = ftell(fid);
        fseek(fid,location,'bof');
        fprintf(fid,'\n'); 
        fseek(fid,-1,'cof'); 

        fprintf(fid,'%s',new_line); 
        fprintf(fid,'\r'); 
        fclose(fid);

        % Remove the job lock
        delete([parent_dir,'runfile_lock']);
        fprintf('Lock file deleted.\n\n');

        previous_data_path = data_path;
        display(['Data image #',int2str(z),' analysis complete'])
    end
end

display('Analysis Complete... Clearing Workspaces')
clear global
clear