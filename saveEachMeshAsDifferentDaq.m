function outputDAQ = saveEachMeshAsDifferentDaq(daq)

if ~exist('daq','var'); daq = []; end
[DAQ,isGlobal] = useLocalOrGlobalDaqFiles(daq);


if ~isfield(daq,'gpopsOutput')
    error('No gpopgOutput field found')
    return
end

output = daq.gpopsOutput;
for i = 1:length(output.meshhistory)
    tempDaq.gpopsOutput.result = daq.gpopsOutput.meshhistory(i).result;
    outputDAQ{i} = parseGpops2toDaq(tempDaq.gpopsOutput,daq.gpopsNames.states,daq.gpopsNames.controls,daq.gpopsNames.indepVar,daq.gpopsNames.units,daq.gpopsNames.names,'solutionField','solution');
end

if isGlobal
    DAQ = [DAQ; outputDAQ];
end