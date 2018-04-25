function DAQ2 = assembleNewDaqAtIndicies(ind,daq,varargin)
%For a given set of indicies, this funtion will assemble a new daq at those
%indicies. This will be useful for comparing specific events.
%INPTUS:
%   ind - either vector (for one daq file) or cell array of indicies for
%       multiple files to make new daq files from class double 2x1 or cell
%       arry of double 2x1, length(ind) == length(DAQ).
%   daq (optional) - pass in a local copy of a daq file or DAQ files
%   varargin - used to set defaults variable
%
%OUTPUTS:
%   DAQ2 - new daq of the files at specific indicies      class struct if
%   one file or claass cell if multiple files
%
%Creation: 21 Sep 2016 - Jeff Anderson
%Future work - at "Indicies of interset" - could make it work if passing in
%"Lap 1" or specific chars.
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% Deal with inputs
if ~exist('daq','var'); daq = []; end
DAQ = useLocalOrGlobalDaqFiles(daq);

%Defaults
defaults = {'newFilename','useOldFilename'                                 %Default is to use same filename as original daq, otherwise pass in new one here
            'newShortFilename','useOldShortFilename'                       %Default is to use same short as original daq, otherwise pass in new one here
            'returnCell',false                                             %Optionally always return a cell, I needed the local daq struct for several other files
            'normalizeIndepVarChannels',true                               %For plotting laps easily
            'channels',[]};                                                %For manually choosing the channels, leave empty for all [] = daqChannels(daq)
setDefaultsForVarargin(defaults,varargin)

%Put ind in a cell if it's a double 2x1
if isa(ind,'double')
    temp = ind;
    clear ind
    ind{1} = temp;
end

%% Operate on files

for iFile = 1:length(DAQ)
    daq = DAQ{iFile};
    newDaq = daq;
    
    if isempty(channels)
        channels = daqChannels(newDaq);
    else %need to remove non-needed channels
        allChannels = daqChannels(newDaq);
        newDaq.rawData = rmfield(newDaq.rawData,setdiff(allChannels,channels));
        
    end
    
    %Indicies of interest - future, could possibly use chars to look for
    %lap 1, etc.
    indicies = ind{iFile};
    
    %Grab just the ind of interest
    for iCh = 1:length(channels)
        channel = channels{iCh};
        try 
            newDaq.rawData.(channel).meas = daq.rawData.(channel).meas(indicies);
            if (strcmp(channel,'time') || strcmp(channel,'distance')) && normalizeIndepVarChannels
                newDaq.rawData.(channel).meas = newDaq.rawData.(channel).meas - newDaq.rawData.(channel).meas(1);
            end
        catch
            warning('Channel "%s" not correct type or dimension for new file on File: %s',channel,displayDaqFiles(daq,'suppressOutput',true));
        end
    end
    
    
    %Strip anything other than header or rawData
    fields = fieldnames(newDaq);
    for iField = 1:length(fields)
        field = fields{iField};
        if strcmp(field,'header') && strcmp(field,'rawData');
            newDaq = rmfield(newDaq,field);
        end
    end
    
    %Deal with filenames
    if ~strcmp(newFilename,'useOldFilename')
        newDaq.header.filename = newFilename;
    end
    
    if ~strcmp(newShortFilename,'useOldShortFilename')
        newDaq.header.shortFilename = newShortFilename;
    end

    
    DAQ2{iFile} = newDaq;
end

if numel(DAQ2) == 1 && ~returnCell
    DAQ2 = newDaq;
end