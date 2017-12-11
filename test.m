%load data from text file
runfile = 'C:\Dropbox\Shared Folders\LLLab Shared\July 2011 Data\Case Data\runfile.txt';
read_fid = fopen(runfile,'r');

z = 1;
current_line = fgetl(read_fid);
while ischar(current_line)
    % Define the new line
    [token,remain] = strtok(current_line,char(9));
    if strcmp(remain,[9,'0'])
        new_line = [token,9,'1'];
        % Opening the file to both read and write 
        fid = fopen(runfile,'r+'); 
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
    end
    
    %Move to the next line
    current_line = fgetl(read_fid);
    z = z + 1;
end    

fclose(read_fid);