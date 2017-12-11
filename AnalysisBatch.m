function AnalysisBatch(parent_dir,data_dir)
%This just serves as a template for analysis scripts.
%Parent_dir is the directory of the data run
%data_dir is the subdirectory of the data you want to analyze e.g.
%'Retrieval' or 'Results'

%% Set variables here
analysis_name = 'ElevHist';
data_type = 'mat';

%% Batch processing structure

if ~exist('analysis_name','var')
    error('Invalid analysis_type, check available functions')
end

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
if ~exist([parent_dir,'Analysis/'],'dir')
    display('Output directory does not exist')
    display('Creating output directory')
    mkdir(parent_dir,'Analysis/');
end

%Check for the subdirectory
if ~exist([parent_dir,'Analysis/',analysis_name,'/'],'dir')
    display('Output sub-directory does not exist')
    display('Creating output sub-directory')
    mkdir(parent_dir,'Analysis/',analysis_name,'/');
end

out_dir = ['Analysis/',analysis_name,'/'];

%Create a list of files in the directory, if necessary
CreateRunfile([parent_dir,data_dir],[parent_dir,out_dir],data_type);

display('Importing File List')
read_fid = fopen([parent_dir,out_dir,'runfile.txt'],'r');

display('Begin batch analysis')
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
        
        % ###### Run the algorithms and save the outputs #####
    
        [output] = ElevHist(full_path);
        save([parent_dir,out_dir,file_name,'.mat'],output)
        
        
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

display('Analysis Complete')