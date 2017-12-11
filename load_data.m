function data_image = load_data(directory,filename)
%% load_data is a function used to import any data image into matlab
% Since data images can be stored in files with many different types of
% file extensions, creating a program flexible enough to handle all
% different types becomes tedious. This script supports:
%
%   - Single Variable Matlab Workspaces (.mat),
%   - Tiff Images (.tiff or .tif)
%
%   ARGUMENTS:
%       directory:
%           -string of the absolute of relative location of data image
%           -e.g. 'C:\Users\Bob\' or 'data\images\'
%       filename:
%           -should be the name of the file you wish to import
%           -should contain the proper file extension in name
%           -e.g. 'image.tif' or 'image.mat'
%           
%   OUTPUTS:
%       data_image
%           -matlab local variable of type imported from file
%
%   AUTHOR: Stevens Institute of Technology - Light and Life Lab
%           Steven Englehardt - englehardt@gmail.com
%
%   VERSION: 
% 1.0 - July 2010
% 1.1 - June 2011
%   -made directory input more flexible

%% Check to see if directory ends with '\'
if directory(length(directory)) ~= '\' && directory(length(directory)) ~= '/'
    directory = [directory,'/'];
end

%% Figure out what the file extension is
remain = filename;
while true
    [type, remain] = strtok(remain, '.'); %#ok
    if isempty(type),  break;  end
    file_type = type;
end   

%% Load import the data, depending on file extension
if strcmpi(file_type,'mat')
    vars = whos('-file', [directory,filename]);
    data_struct = load([directory,filename]);
    data_image = data_struct.(vars(1).name);
elseif strcmpi(file_type,'tif') || strcmpi(file_type,'tiff')
    data_image = imread([directory,filename],'tif');
elseif strcmpi(file_type,'bmp')
    data_image = imread([directory,filename],'bmp');
else
    error('This type of file extension is not supported');
end
   