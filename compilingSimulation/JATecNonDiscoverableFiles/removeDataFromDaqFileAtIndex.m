function daq = removeDataFromDaqFileAtIndex(indicies,daq,varargin)
%This function will remove data from the daq file passed in at the indicies
%ind that are passed in.
%INPUTS:
%    indicies - vector of indicies to remve from the daq file              class int nx1
%          OR a cell array of vectors for each file                        class cell (length(daq)x1)
%    daq - local instance of the daq file to remove indicies from          class struct
%OUTPUT:
%    daq - final version of the daq file                                   class struct
%
%Creation: 22 Sep 2015 - Jeff Anderson
%
%Updates: 12 Sept 2016 - Matthew Pearson and Jeff Anderson
%                           Applies global DAQ is no daq is passed in
%Updates: 19 Jun 2017 - Jeff Anderson - I made the notes optionall (not by
%    default) add a note of each index removed.
%Updated:  02 Aug 2017 - Jeff Anderson, delt with removing indicies from
%    lapInfo section.
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

defaults = {'addNotesAboutWhichIndWereRemoved' false                       %This is sometimes slow and I rarely need to know what original indicies were removed
            'noteAboutWhyDataWasRemoved','data was removed from this channel'
            'saveOrigMeasureIfNotPresent',true};
setDefaultsForVarargin(defaults,varargin)

if ~exist('daq','var'); daq = []; end
DAQ = useLocalOrGlobalDaqFiles(daq);




for iFile = 1:length(DAQ)
    if isa(indicies,'cell')
        ind = indicies{iFile};
    else
        ind = indicies;
    end

    %Sort them
    ind = sort(ind,'descend');
    
    daq = DAQ{iFile};
    channels = daqChannels(daq);
    
    
    %Deal with lapInfo
    if isfield(daq,'lapInfo')
        ch = getIndependantVariableChannel(daq);
        getChannelDataFromDaqFile(daq,{'s' ch})
        lapStatus = zeros(size(s));
        for iLap = 1:length(daq.lapInfo.laps)
            lapInd = [daq.lapInfo.laps{iLap}.ind(1):daq.lapInfo.laps{iLap}.ind(end)];
            lapStatus(lapInd) = iLap;
        end
        
        %Now remove all ind that need to be removed and recalculate the lap
        %indicies
        lapStatus(ind) = [];
        for iLap = 1:length(daq.lapInfo.laps)
            lapInd = find(lapStatus == iLap);
            daq.lapInfo.laps{iLap}.ind = [min(lapInd) max(lapInd)];
        end
    end
    
    %Loop through the channels
    for iCh = 1:length(channels)
        %If there is no backup of the original data
        if ~isfield(daq.rawData.(channels{iCh}),'originalMeas') && saveOrigMeasureIfNotPresent
            daq.rawData.(channels{iCh}).originalMeas = daq.rawData.(channels{iCh}).meas;
        end

        daq.rawData.(channels{iCh}).meas(ind) = [];

        %Add some notes if none exist
        if addNotesAboutWhichIndWereRemoved
            notes = 'Indicies';                                                    %Start the comments by saying which indices we are removing
            for i = 1:length(ind)
                notes = [notes sprintf(' %i,',ind(i))];                            %#ok dyn growth...which ones
            end
            notes(end) = [];                                                       %remove last comma;
            notes = [notes sprintf(' were removed on %s',datestr(now,'YYYY-mm-dd'))];           %#ok dyn growth...add date                          
            
        else
            notes = noteAboutWhyDataWasRemoved;
        end
        daq.rawData.(channels{iCh}) = addNotesToDaqFile(daq.rawData.(channels{iCh}),notes); %save it
    end
    

    DAQ{iFile} = daq;
end    