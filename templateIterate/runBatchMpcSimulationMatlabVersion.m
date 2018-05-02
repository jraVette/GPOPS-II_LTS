
%Add the matlab code to the path
[templateDirectoryPath,~,~] = fileparts(mfilename('fullpath'));            %Current directory
indFileSep = strfind(templateDirectoryPath,filesep());                     %Remove last directory, grab the file seperators
simRootDirectory = templateDirectoryPath(1:indFileSep(end-1)-1);             %Get the address to the last filesep
codeRepositoryPath = fullfile(simRootDirectory,'compilingSimulation');     %Now add the 'compilingSimulation' directory
addpath(genpath(codeRepositoryPath));

%Load the instance and run
load('daqFile.mat')
gpopsMPC(daq)

