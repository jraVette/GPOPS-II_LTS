function [allChannels, allCommonChannels, unCommonChannelsPerFile, localChannels, allUnits, allCommonUnits, unCommonUnitsPerFile, localUnits] = daqChannels(daq)
% [allChannels, allCommonChannels, unCommonChannelsPerFile, localChannels, allUnits, allCommonUnits, unCommonUnitsPerFile, localUnits] = daqChannels(daq)
%This funciton will return all common channels that are loaded into the DAQ
%structure. 
%
%INPUTS:
%    daq (OPTIONAL) - either a single instance or cell arrary of daq files
%                     to look at.
%
%OUTPUTS:
%    allChannels = list of all channels available            class row cell 
%    allCommonChannels = list of all channels available on every file
%                                                            class row cell
%    unCommonChannelsPerFile - cell vector (size nFilesx1) with a cell
%                              vector in each discribing the channels that
%                              aren't common.                class row cell
%    localChannels - are a channel list for each file        class row cell
%    units - repeat of the above 4 outputs but units
%
%EXAMPLE:
%    channels = daqChannels(); will return all common channels
%    channels = daqChannels(DAQ{1}); will return the channels on the file
%                                    DAQ{1}.
%
%Creation: 28 May 2013 - Jeff Anderson
%Updated:  05 Feb 2015 - Jeff Anderson - edited to pass in a specific file
%                                        or use all laoded files
%Updated:  14 May 2015 - Jeff Anderson - fixed to show all common channels,
%                                        all channels, and the uncommon 
%                                        ones per file
%Updated:  27 Aug 2015 - Jeff Anderson - Fixed to error out appropiately if
%                                        no daq file is loaded
%Updated   10 Dec 2015 - Jeff Anderson - Added support for units
%Updated:  13 Feb 2017 - Jeff Anderson
%    I only checked the common/uncommon channels if the daq passed was
%    multiple daq instances. Using ismember can be slow
%Updated:  10 Jun 2017 - Jeff Anderson - bug fixes for only one file passed
%    in.
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% Deal with inputs
if ~exist('daq','var')
    daq = [];
end
[allDaqFiles,~] = useLocalOrGlobalDaqFiles(daq);

%If only one file, we don't need to find differences
if length(allDaqFiles)>1
    checkUncommon = true;
else
    checkUncommon = false;
end


%% Deal w/ each file and get the local channels and units
allChannels = [];
allUnits = [];
localChannels = cell(length(allDaqFiles),1);
localUnits = cell(length(allDaqFiles),1);
for iFile = 1:length(allDaqFiles)                                          %Loop through files
    daq = allDaqFiles{iFile};                                              %Grab local instance
    if ~isfield(daq,'rawData'); 
        channelsOnThisFile = [];
        warning('No rawData found on file: %s',displayDaqFiles(daq,'suppressOutput',true))
    else
        channelsOnThisFile = fieldnames(daq.rawData);                      %Grab local channel names
        for iCh = 1:length(channelsOnThisFile)                             %Loop through the channels
            if isfield(daq.rawData.(channelsOnThisFile{iCh}),'units')      %Try to get the units
                unitsOnThisFile{iCh} = daq.rawData.(channelsOnThisFile{iCh}).units;
            else
                unitsOnThisFile{iCh} = '';
            end
        end
    end
    unitsOnThisFile = cell(length(channelsOnThisFile),1);                  %Make a cell array for the units

    
    allChannels = [allChannels' channelsOnThisFile']';                     %Add this to the full list                   
    allUnits    = [allUnits'    unitsOnThisFile']';                        %Add units to the full list
    
    localChannels{iFile} = channelsOnThisFile;                             %Save the local intance
    localUnits{iFile} = unitsOnThisFile;                                   %Save the local units
end

%% Make the unique list of all channels available
[allChannels,ind,~] = unique(allChannels);                                 %What are the channels available
allUnits = allUnits(ind);


%% Get all common channels and units
if checkUncommon
    flagOnAllChannels = ones(length(allChannels),1);                           %First assume it's on all the files
    for iCh = 1:length(allChannels)                                            %Loop through channles
       channelToCheck = allChannels{iCh};                                      %Grab the channel to check all files for
       for iFile = 1:length(localChannels)                                     %Loop through the files
           if ~ismember(channelToCheck,localChannels{iFile})                   %If that channel is not on the file...
               flagOnAllChannels(iCh) = 0;                                     %Flag it     
           end
       end                                                                     %End looping through files
    end                                                                        %End looping through channles
    ind = find(flagOnAllChannels);                                             %Get teh index of the common channels
    allCommonChannels = cell(length(ind),1);                                   %Preallocate the common channels
    allCommonUnits    =  cell(length(ind),1);
    for i = 1:length(ind)                                                      %Loop through them and save
        allCommonChannels{i} = allChannels{ind(i)};
        allCommonUnits{i}    = allUnits{ind(i)};
    end
    

    %% Now see what are the uncommon channels per file loaded
    unCommonChannelsPerFile = cell(size(allDaqFiles));
    unCommonUnitsPerFile    = cell(size(allDaqFiles));
    if exist('localChannels','var')
        for iFile = 1:length(localChannels)
           [unCommonChannelsPerFile{iFile},ind] = setdiff(localChannels{iFile},allCommonChannels); 
           unCommonUnitsPerFile{iFile} = localUnits{iFile}(ind);
        end
    else
        error('No Channels Found, please load DAQ files')
    end
else
    allCommonChannels = allChannels;
    unCommonChannelsPerFile = cell(length(allDaqFiles));
    allCommonUnits = allUnits;
    unCommonUnitsPerFile = cell(length(allDaqFiles));
    
    
end

