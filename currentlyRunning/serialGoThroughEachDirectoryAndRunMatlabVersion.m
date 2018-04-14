directories = dir('*DOE*');
currentDirectory = pwd;
for iDir = 1:length(directories)
    cd(directories(iDir).name)
    runBatchMpcSimulationMatlabVersion
    cd(currentDirectory)
end