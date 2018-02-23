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

%Load up the last one to see if we have an answer
for iHorizon = max(horizonNumber):-1:min(horizonNumber)
    daqToLoad = load(files(iHorizon).name);
    if ~isempty(daqToLoad.segDaq)
        DAQ{1} = daqToLoad.masterDaq;
        break
    else
        files(iHorizon) = [];
        horizonNumber(iHorizon) = [];
    end
end


for iFile = 2:length(horizonNumber)+1
    daqToLoad = load(files(iFile-1).name);
    fields = fieldnames(daqToLoad);
    DAQ{iFile} = daqToLoad.segDaq;
    
    if fixDistance
        currentDistance = daqToLoad.masterDaq.status.currentDistance;
        DAQ{iFile}.rawData.distance.meas = DAQ{iFile}.rawData.distance.meas+currentDistance;
    end
    
    
    for iField = 1:length(fields)
        field = fields{iField};
        if ~strcmp(field,'segDaq')
            DAQ{iFile}.(field) = daqToLoad.(field);
        end
    end
    
    DAQ{iFile}.header.shortFilename = sprintf('Horizon%03i',horizonNumber(iFile-1));
    
end

