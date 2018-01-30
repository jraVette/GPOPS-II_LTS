global DAQ
files = dir('*Horizon*.mat')
for iFile = 1:length({files.name}')
    load(files(iFile).name);
    DAQ{iFile} = segDaq;
end