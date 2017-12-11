function CreateRunfile(parent_dir,ext,start_index,end_index)
% CreateRunfile creates a text file containing lists of all images to
% include in the data processing, a binary value indicating if the files
% has been run or not, and a run time if run.

%parent_dir: The directory where data is located e.g. 'C:/data/case0'
%ext: extension for files to be analyzed e.g. 'tif'
%start_index (end_index): image index number to start (end) processing at
%   e.g. to process images 100 - 200 use CreateRunfile(...,...,100,200)
%   TO PROCESS ALL IMAGES USE: CreateRunfile(...,...,1,-1)

% Check to see if directory ends with a '\'
[part,remainder] = strtok(parent_dir,'/\');
if isunix
    parent_dir = ['/',part];
else
    parent_dir = part;
end

if isunix
    while ~isempty(part)
        [part,remainder] = strtok(remainder,'/\');
        parent_dir = [parent_dir,'/',part];
    end
else
    while ~isempty(part)
        [part,remainder] = strtok(remainder,'/\');
        parent_dir = [parent_dir,'\',part];
    end
end

%Get a recursive file list
display('Retrieving list of all images...')
if isunix
    File_List = rdir([parent_dir,'**/*.',ext]);
else
    File_List = rdir([parent_dir,'**\*.',ext]);
end

%Open the output file
if exist([parent_dir,'runfile.txt'],'file')
    display('Runfile already exists, continuing from there.')
    return;
else
    fid = fopen([parent_dir,'runfile.txt'],'w');
end

previous_data_path = '';
index = 1;

%Run through the file list, saving the name and path of each file, as well
%as a bool false value
for z = 1:length(File_List)
    
    %Parse the full path down to just the file name
    if isunix
        data_path = '/';
    else
        data_path = '';
    end
    token = '';
    remain = File_List(z).name;
    while ~isempty(remain)
        if ~isempty(token)
            data_path = [data_path,token,'/'];
        end
        [token,remain] = strtok(remain,'/\');
    end
    
    if (~strcmp(data_path,previous_data_path))
        index = 1;
        display(['Listing images from dir: ',data_path])
    end
    
    %Save if we are in the right range
    if( (end_index == -1 && index >= start_index) || (end_index >= index && start_index <= index) )
        token = '';
        remain = File_List(z).name;
        if isunix
            fprintf(fid,'/');
        end
        while ~isempty(remain)
            if ~isempty(token)
                fprintf(fid,[token,'/']);
            end
            [token,remain] = strtok(remain,'/\');
        end
        fprintf(fid,token);
        fprintf(fid,'\t0\r');
    end
    
    previous_data_path = data_path;
    index = index + 1;
    
end

%Close the output file
fclose(fid);