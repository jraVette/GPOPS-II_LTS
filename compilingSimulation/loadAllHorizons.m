%Script to load all horizons into a daq
%Created 02 Feb 2018 - Jeff Anderson
%Updated 27 Feb 2018 - Jeff Anderson - updated to new method
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
c
global DAQ

fixDistance=true;
files = dir('*Horizon*.mat');

currentColorOrder = plotColor;
plotColor([],[[0 0 0];currentColorOrder;currentColorOrder]);
currentColorOrder = plotColor;
%First find all horizon numbers
for iFile = 1:length(files)
    horizonKeyWord = 'Horizon';
    indHorizon = strfind(files(iFile).name,horizonKeyWord);
    horizonNumber(iFile) = str2num(files(iFile).name(indHorizon+length(horizonKeyWord):indHorizon+length(horizonKeyWord)+2));
end
horizonNumber = unique(horizonNumber);

%Load the last master daq possible
for iHorizon = max(horizonNumber):-1:min(horizonNumber)
    filename = sprintf('Horizon%03i-MasterDaq.mat',horizonNumber(iHorizon));
    if exist(filename,'file')
        load(filename)
        fprintf('Loading Master Daq: %s\n',displayDaqFiles(daq,'suppressOutput',true,'lookForShortFilename',false));
        DAQ{1} = daq;
        break
    end
end

%Load all the horizons
for iHorizon = 1:length(horizonNumber)
    fprintf('Loading horizon %03i\n',horizonNumber(iHorizon));
    filename = sprintf('Horizon%03i-OCP.mat',horizonNumber(iHorizon));
    load(filename)
    
    if fixDistance && isfield(daq,'rawData')
        daq.rawData.distance.meas = daq.rawData.distance.meas+daq.gpopsSetup.auxdata.currentDistance;
    end
    
    daq.header.shortFilename = sprintf('Horizon%03i',horizonNumber(iHorizon));
    DAQ{end+1} = daq;
    
end

