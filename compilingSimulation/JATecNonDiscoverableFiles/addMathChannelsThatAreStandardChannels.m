function rawDataStruct = addMathChannelsThatAreStandardChannels(rawDataStruct,standardChannelName,channelData,notes,units,varargin)
%This function will try to use the standardNames.xls to port in the
%stanadard information for the channel so user won't have to type it in
%again. Uses a single channel at a time!!!
%INPUTS:
%   rawDataStruct - structure of all the raw data (place to add channels)
%                                                              class struct
%   standardChannelName - standard channel name, must be found in
%                         channelNames.xls, col C              class char
%   channelData - row vector of data to add                    class double
%   notes (optional) - cell array of notes to add to channel   class char
%   units (optional) to overide units in channelNames.xls      class char
%OUTPUS:
%   rawDataStruct - updated rawData structure                  class struct
%
%Example:
%   >> daq.rawData = addMathChannelsThatAreStandardChannels(daq.rawData,'time',0:0.1:1);
%   This will add a time channel and pull in info from
%   channelInfo.xls
%
%Creation: 12 Apr 2016 - Jeff Anderson
%Updated:  01 Oct 2017 - Jeff Anderson - option to suppress warnings
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

defaults = {'suppressWarning' false
            'suppressOutput'  false};
setDefaultsForVarargin(defaults,varargin);

%Inputs
if ~exist('notes','var'); notes = ''; end
if ~exist('units','var'); units = ''; end %This is pulled 

%Get the channel names
[~,commonHeaders]    = parseChannelNamesIntoMat();

if ~ismember(standardChannelName,fieldnames(commonHeaders))                        %if the common channel isn't found
    if ~suppressWarning
        warning('Channel %s, not a standard name!',standardChannelName);
    end
    rawDataStruct.(standardChannelName) = createDaqChannelData(channelData,units,standardChannelName,'mathChannel',true);
else                                                                   %Ok, we found the channel, port it in
    commonChannelInfo = commonHeaders.(standardChannelName);
    fields = fieldnames(commonChannelInfo);

    for iField = 1:length(fields)                                      %Loop through and port the fields over, but don't port aliases, dont' need it
        if ~strcmp(fields{iField},'aliases')
            rawDataStruct.(standardChannelName).(fields{iField}) = commonChannelInfo.(fields{iField});
        end
    end

    %Port notes in
    rawDataStruct.(standardChannelName) = addNotesToDaqFile(rawDataStruct.(standardChannelName) ,notes);

    %If we want to override units
    if ~strcmp(units,'')
        rawDataStruct.(standardChannelName).units = units;
    end

    %consider this a math channel
    rawDataStruct.(standardChannelName).mathChannel = true;

    %Add data
    rawDataStruct.(standardChannelName).meas = rowVector(channelData);
    if ~suppressOutput
        fprintf('Channel: %s added to daq file\n',standardChannelName);
    end
    
    
end





