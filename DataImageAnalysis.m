function DataImageAnalysis(parent_dir)
%This script provides the following statistics on full sets of
%images
%FOR EACH CHANNEL (Color)
%    Histogram
%    Mean
%    Variance
%    Maximum
%    Minimum

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

%Check to see if the root contains images, and if any subdirs do
dir_index = 1;
File_List = dir([parent_dir,'/*.tif']);
if ~isempty(File_List)
    directory_list{dir_index} = parent_dir;
    dir_index = dir_index + 1;
end

File_List = dir(parent_dir);
for z = 3:length(File_List)
    if File_List(z).isdir
        New_File_List = dir([parent_dir,File_List(z).name,'/*.tif']);
        if ~isempty(New_File_List)
            %has tif images, add to list of dirs
            directory_list{dir_index} = [parent_dir,File_List(z).name,'/'];
            dir_index = dir_index + 1;
        end
    end
end

%Check for output directory existance
if ~exist([parent_dir,'Analysis/'],'dir')
    display('Output directory does not exist')
    display('Creating output directory')
    mkdir(parent_dir,'Analysis/');
end

clear File_List
for z = 1:size(directory_list,2)
    
    display(['Processing in directory: ',directory_list{z}])
    
    %Initialize storage variables
    total = zeros(1024,1280,3);
    total_square = zeros(1024,1280,3);
    total_mean = [0;0;0];
    total_var = [0;0;0];
    total_hist = zeros(3,1024);
    maximum = [-Inf;-Inf;-Inf];
    minimum = [Inf;Inf;Inf];
    
    %Check for output directory existance
    sub_dir = strrep(directory_list{z},parent_dir,'');
    if ~exist([parent_dir,'Analysis/',sub_dir],'dir')
        display('Output directory does not exist')
        display('Creating output directory')
        mkdir([parent_dir,'Analysis/',sub_dir]);
    end
    
    %Load correction matrices
    %display('Loading new gain.mat and offset.mat')
    %load([directory_list{z},'/gain.mat']);
    %load([directory_list{z},'/offset.mat']);
    
    File_List = dir([directory_list{z},'/*.tif']);
    
    for k = 1:length(File_List)
    
        %Load the image - mean, max, min
        data = double(demosaic(imread([directory_list{z},File_List(z).name],'tif'),'gbrg'));
        %data = double(demosaic(uint16((double(data) - offset)),'gbrg')).*gain;
        
        k
        
        for i = 1:3
            total(:,:,i) = total(:,:,i) + data(:,:,i);
            total_square(:,:,i) = total_square(:,:,i) + data(:,:,i).^2;
            maximum(i) = max(maximum(i),max(max(data(:,:,i))));
            minimum(i) = min(minimum(i),min(min(data(:,:,i))));
            total_hist(i,:) = total_hist(i,:) + hist(reshape(data(:,:,i),1024*1280,1),0:1:1023);
        end
    end
    
    for i = 1:3
        total_mean(i) = mean(mean(total(:,:,i)))/length(File_List);
        total_var(i) = (1/length(File_List)) * mean(mean(total_square(:,:,i))) - total_mean(i)^2;
    end
    
    %Save the results in a mat file
    save([parent_dir,'Analysis/',sub_dir,'analysis_results.mat'],'total_mean','total_var','maximum','minimum');
    figure; bar(0:1:1023,total_hist(1,:)')
    title('Red Histogram')
    saveas(gcf,[parent_dir,'Analysis/',sub_dir,'red_hist.png'])
    close(gcf)
    figure; bar(0:1:1023,total_hist(2,:)')
    title('Green Histogram')
    saveas(gcf,[parent_dir,'Analysis/',sub_dir,'green_hist.png'])
    close(gcf)
    figure; bar(0:1:1023,total_hist(3,:)')
    title('Blue Histogram')
    saveas(gcf,[parent_dir,'Analysis/',sub_dir,'blue_hist.png'])
    close(gcf)
    
end