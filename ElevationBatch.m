function ElevationBatch(parent_dir)
%% ElevationBatch runs the ElevationRetrieval routine on a directory of data
% The data is assumed to be in tif form, and will be output to a
% subdirectory \Output\ with .mat files of the same names as the input
% images. The mat files contain the final surface normals and surface
% elevation matrices.

%Lock-File implementation inspired by: http://www.mathworks.com/matlabcentral/newsreader/view_thread/287163

% Generate a process ID
c = clock;
myID = num2str(etime(c,[2012 0 0 0 0 0]));
complete = false;

%% Load and prepare flat surface
display('Loading constants from data directory')

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



%Import the list of files from that directory
display('Importing File List')
read_fid = fopen([parent_dir,'runfile.txt'],'r');
    
%Load the model and system constants
load([parent_dir,'model_constants.mat'])
system_constants;
tic;
disp('beginning prefactorization')
[A R W e] = R_factor(1024,1280);
disp(['factorization time ' num2str(toc)]);
iteration_vector = [2,1];

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
z = 1;
current_line = fgetl(read_fid); %Load first line of data

while ischar(current_line) %Reads until end of file
   
    % Define the new line
    [full_path,flag] = strtok(current_line,char(9));
    
    if strcmp(flag,[9,'0']) %We haven't analyzed this data yet
        
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
        
        %Load the image
        data = load_data(data_path,file_name);
        data = double(demosaic(uint16((double(data) - offset)),'gbrg')).*gain;
        
        display([10,'Data Image: ',file_name])
        display(['Location: ',data_path,10])
        
        %Iterate through retrieval
        [surface] = SolveSystem(data,u_fp,X,Y,PL,spectra_cube, ...
            flat_surf,fp_res,iteration_vector,z_fp,W,A,R,e); %#ok
        surface(:,:,3) = surface(:,:,3) - mean(reshape(surface(:,:,3),...
            numel(surface(:,:,3)),1));
    
        elevation = single(surface(:,:,3));
    
        %Check for output directory existance
        if ~exist([parent_dir,'Retrieval/',sub_dir],'dir')
            display('Output directory does not exist')
            display('Creating output directory')
            mkdir([parent_dir,'Retrieval/',sub_dir]);
        end
    
        %Save x-y positions if this is the first run
        if z == 1
            xy_position = surface(:,:,1:2);
            save([parent_dir,'Retrieval/',sub_dir,'xy_position.mat'],'xy_position')
        end
    
        %Save the mat file
        save([parent_dir,'Retrieval/',sub_dir,strtok(file_name,'.'),'.mat'],'elevation')
        
        %Update the run file
        new_line = [data_path,file_name,9,'1'];
        % Opening the file to both read and write 
        fid = fopen([parent_dir,'runfile.txt'],'r+'); 
        % Can change 'loc' to insert data at any line in file. 
        % When loc=2, data will be inserted at line 3 
        loc = z-1; 
        for i = 1:loc 
            %Used FGETL to move file pointer a whole line at a time: 
            %See FGETL section below for more information 
            temp_line = fgetl(fid);
        end; 
        location = ftell(fid);
        fseek(fid,location,'bof');
        fprintf(fid,'\n'); 
        fseek(fid,-1,'cof'); 

        % This call to FPRINTF utilizes the vectorized feature of the 
        % MATLAB version: see below for more information 
        fprintf(fid,'%s',new_line); 
        fprintf(fid,'\r'); 
        fclose(fid);

        display(['Data image #',int2str(z),' analysis complete'])
    end
    
    %Move to the next line in runfile
    current_line = fgetl(read_fid);
    z = z + 1;
end
