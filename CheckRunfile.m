function CheckRunfile(parent_dir)
% CheckRunfile checks the runfile to make sure all data
% was processed successfully. If it finds any instances as failed
% it marks them as unprocessed and reports this to the user.

%parent_dir: The directory where data is located e.g. 'C:/data/case0'

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

%Open the output file
if ~exist([parent_dir,'runfile.txt'],'file')
    display('Runfile does not exist! Unable to determine state of processing.')
    return;
else
    read_fid = fopen([parent_dir,'runfile.txt'],'r');
end

%Run through the runfile
unprocessed_count = 0;
failed_count = 0;
successful_count = 0;
z = 1;
current_line = fgetl(read_fid); %Load first line of data
while ischar(current_line) %Reads until end of file
   
    % Define the new line
    [full_path,flag] = strtok(current_line,char(9));
    
    if strcmp(flag,[9,'0']) %We haven't analyzed this data yet
        unprocessed_count = unprocessed_count + 1;
    elseif strcmp(flag,[9,'2']) %Processing failed, correcting
        failed_count = failed_count + 1;
        
        % Update the run file
        new_line = [full_path,9,'0'];
        % Opening the file to both read and write 
        fid = fopen([parent_dir,'runfile.txt'],'r+'); 
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
    else %Successfully Processed Data
        successful_count = successful_count + 1;
    end
    
    %Move to the next line in runfile
    z = z + 1;
    current_line = fgetl(read_fid);
    
end

display([10,'Out of ',num2str(successful_count+unprocessed_count+failed_count),...
    ' data images:', 10, 10, num2str(successful_count), ' have been processed successfully.'])

if failed_count > 0 || unprocessed_count > 0 
    display([num2str(unprocessed_count), ' images were never processed.'])
    display([num2str(failed_count), ' images failed to process successfully, the flags have been reset to 0.'])
    display([10, 'Processing NOT COMPLETE!',10,'Re-run ElevationBatch_global and check again once that has completed'])
else
    display('All images have been processed successfully, processing COMPLETE!')
end

%Close the output file
fclose(read_fid);