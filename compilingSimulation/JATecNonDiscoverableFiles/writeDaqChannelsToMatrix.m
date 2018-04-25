function outputs = writeDaqChannelsToMatrix(daq,varargin)
%This function will output all the channel data to a matrix for
%mainpulation later.
%INPUTS: 
%    (optional) daq - local instance of a daq file             class struct
%
%OUTPUS:
%    outputs - matrix of the output of each file, if daq is used or only
%    one file is loaded, then class(output)=double, otherwise,
%    class(output)=cell where size(output) = cell(length(DAQ),1)
%
%Creation: 1 Oct 2015 - Jeff Anderson
%Updated:  11 Nov 2017 - Jeff Anderson - added just select channels
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
selectedChannels = [];
%Deal with inputs
defaults = {'selectedChannels', [] %Cell array of channels
            'wheelPositionSuffix',[]}; %Wheel position to loop over
setDefaultsForVarargin(defaults,varargin);

if ~exist('daq','var');    daq = [];  end
DAQ = useLocalOrGlobalDaqFiles(daq);

%Make sure selected channels are a cell
if ~isempty(selectedChannels) && ~isa(selectedChannels,'cell')
    temp = selectedChannels; clear selectedChannels
    selectedChannels{1} = temp;
end

%Wheel positions
count = 1;
if ~isempty(wheelPositionSuffix) %#ok defined in defaults
    %Make sure wheel positions are a cell
    if isa(wheelPositionSuffix,'char'); temp = wheelPositionSuffix; clear wheelPositionSuffix; wheelPositionSuffix{1} = temp; end
    
    %Loop through selected channels and make a new list
    for iCh = 1:length(selectedChannels)
        for iPos = 1:length(wheelPositionSuffix)
            newSelectedChannels{count} = [selectedChannels{iCh} wheelPositionSuffix{iPos}]; %#ok dyn growth
            count = count+1;
        end
    end
    selectedChannels = newSelectedChannels;
end


%Grab the channel data
for iFile = 1:length(DAQ)
    daq = DAQ{iFile};
    if isempty(selectedChannels)
        selectedChannels = daqChannels(daq);
    end
    
    [nRow,nCol] = size(daq.rawData.(selectedChannels{1}).meas);
    if nCol ~= 1
        error('Channel %s on file %s is not a row vector, please correct',indepCh,daq.header.filename);
    end
    tempMatrix = zeros(nRow,length(selectedChannels));
    
    
    for iCh = 1:length(selectedChannels)
        tempChData = rowVector(daq.rawData.(selectedChannels{iCh}).meas);
        [r,c] = size(tempChData);
        if r ~= nRow || c ~= nCol %If it isn't the same as the indp var, then error out
            %Firt try to transpose
            error('Channel %s on file %s is not the same dimenision as %s, please correct. Cannot continue',selectedChannels{iCh},daq.header.filename,indepCh);
        end
        
        tempMatrix(:,iCh) = tempChData;
    end
    
    outputs{iFile} = tempMatrix;
end

%make it not a cell if it is the only file
if length(outputs) == 1
    outputs = outputs{1};
end
