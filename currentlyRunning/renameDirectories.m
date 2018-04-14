c
directories = dir('*DOE*');
currentDirectory = pwd;
for iDir = 1:length(directories)
    cd(directories(iDir).name)
    DAQ = [];
    loadFiles('2018*')
    times = calculateManeuveringTime;
    cd ..
    

    newName = [directories(iDir).name sprintf('_maxIter%i_maxCol%i',daq.header.setup.mesh.maxiterations,daq.header.setup.mesh.colpointsmax)];
    newName = [newName '_' num2str(times)];
    newName = strrep(newName,'.',',');
    sysCommand = sprintf('mv %s %s',directories(iDir).name,newName);
    system(sysCommand)
    
    cd(currentDirectory)
end