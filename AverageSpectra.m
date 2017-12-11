%AverageSpectra.m is used to prepare the spectra data for SpectrumCalc.m

% DIRECTORY and DATA STORAGE
% The spectrum data should be organized into single directory specified by 
% the 'directory' input variable. The 'directory' should have subdirectories
% labeled R#C# corresponding to the row and column of data collection.
% These subdirs should then contain several .dat files with copies of
% several data runs from that same location. Run 'AverageSpectra.m' on the
% parent directory to average these .dat files and convert to .mat files.

if directory(length(directory)) ~= '\'
    directory = [directory,'\'];
end

Dir_List = dir(directory);

% z = 1 will be '.' and z = 2 will be '..'
for z = 3:length(Dir_List)
    
    %Check if a directory
    if Dir_List(z).isdir
        
        %Read in the data files
        File_List = dir([directory,Dir_List(z).name,'\*.txt']);
        clear average
        
        A = dlmread([directory,Dir_List(z).name,'\',File_List(3).name],'\t',[17 0 3664 1]);
        average = A;
        n = 1;
        
        for i = 4:length(File_List)
            A = dlmread([directory,Dir_List(z).name,'\',File_List(i).name],'\t',[17 0 3664 1]);
            average = ((average * n) + A)/(n+1);
            n = n + 1;
        end
        
        eval([Dir_List(z).name,'Spectra = average;'])
        save([directory,Dir_List(z).name,'\',Dir_List(z).name,'Spectra.mat'],[Dir_List(z).name,'Spectra']);
        eval(['clear ',[Dir_List(z).name,'Spectra']])
    end
end