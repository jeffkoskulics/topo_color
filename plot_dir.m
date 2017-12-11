function plot_dir(parent_dir,type)
% This function allows you to specify an output directory and plot the
% images of each elevation

%directory -- should be the directory of .mat outputs that contain the
%             surfaces
%type -- should be either 'fig' for matlab figures or 'tif' for tiff
%        images

if ~strcmp(type,'fig') && ~strcmp(type,'tif')
    error('set type equal to ''fig'' for matlab figures or ''tif'' for images')
end

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

%Check for output directory existance
if isunix
    if ~exist([parent_dir,'Plots/'],'dir')
        display('Output directory does not exist')
        display('Creating output directory')
        mkdir(parent_dir,'Plots/');
    end
else
    if ~exist([parent_dir,'Plots\'],'dir')
        display('Output directory does not exist')
        display('Creating output directory')
        mkdir(parent_dir,'Plots\');
    end
end

%Get a recursive file list
if isunix
    File_List = rdir([parent_dir,'**/image*.mat']);
else
    File_List = rdir([parent_dir,'**\image*.mat']);
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

previous_data_path = '';

for z = 1:length(File_List)
  
    %Parse the full path down to just the file name
    if isunix
        data_path = '/';
    else
        data_path = '';
    end
    token = '';
    file_name = '';
    remain = File_List(z).name;
    while ~isempty(remain)
        if ~isempty(token)
            data_path = [data_path,token,'/'];
            file_name = remain(2:length(remain));
        end
        [token,remain] = strtok(remain,'/\');
    end
    sub_dir = strrep(data_path,parent_dir,'');
    [save_name, remain] = strtok(file_name, '.'); %#ok

    %Check for output directory existance
    if ~exist([parent_dir,'Plots/',sub_dir],'dir')
        display('Output directory does not exist')
        display('Creating output directory')
        mkdir([parent_dir,'Plots/',sub_dir]);
    end
   
    %Load the image
    if (~strcmp(data_path,previous_data_path))
        load([data_path,'xy_position.mat']);
    end
    
    load([data_path,file_name],'elevation');
    plot_surface(double(cat(3,xy_position,elevation)));
    
    title([save_name,' Elevation Plot'])
    
    if strcmp(type,'fig')
        saveas(gcf,[parent_dir,'Plots/',sub_dir,save_name,'_ElevationPlot.fig'],'fig');
    end
    
    if strcmp(type,'tif')
        saveas(gcf,[parent_dir,'Plots/',sub_dir,save_name,'_ElevationPlot.tif'],'tiff');
    end
    
    close(gcf);
    clear elevation
    
    previous_data_path = data_path;
end