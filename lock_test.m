%Lock-File implementation inspired by: http://www.mathworks.com/matlabcentral/newsreader/view_thread/287163

% Generate a process ID
c = clock;
myID = num2str(etime(c,[2012 0 0 0 0 0]));
complete = false;

parent_dir = 'C:\Users\Steven\Desktop\';

while ~complete
    %% Pull the next line of fresh data
    
    % Lock the job list
    lockSuccess = 0;
    while ~lockSuccess
        % Wait for lock-file to disappear
        fprintf('Checking job-list availability...');
        while exist([parent_dir 'runfile_lock'],'file')
            file_info = dir([parent_dir 'runfile_lock']);
            if (~isempty(file_info) && (now - file_info.datenum > .0035) ) %5 mins ~ .0035 in datenum format
                fprintf('Lock file more than 5 minutes old, a crash is likely... Deleting lock file.\n');
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

    %%  Process data
    
    if ~complete
        pause(randi(10))
    
        % Lock the job list
        lockSuccess = 0;
        while ~lockSuccess
            % Wait for lock-file to disappear
            fprintf('Checking job-list availability...');
            while exist([parent_dir 'runfile_lock'],'file')
                pause(0.1);
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
    
    end
end

fprintf('Analysis Complete... Exiting.\n\n');