classdef jatecDaqDataSelector
%This class was written as a backend to getChannelDataFromDaqFile().  See
%creator method for detailed description and usage and use function
%getChannelDataFromDaqFile() to extract daq data.
%
%Written: 6 Jan 2016 - Jeff Anderson
%Updated: 11 Sep 2017 - Jeff Anderson - added special case of getting an
%    index channel from a daq file.
%Updated: 10 Oct 2017 - Jeff Anderson - added option to return origMeas
%    from the daq file.
%Updated: 10 Oct 2017 - Jeff Anderson - origMeas changed to originalMeas per documentation
%Updated: 19 Oct 2017 - Jeff Anderson - bug fixes on unit conversions
%Updated: 23 Feb 2018 - Jeff Anderson. Bug fix, for the case of an invalid
%    file, it won't error out, gracefullly exists now.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    properties 
        daq
        channelsToGetData
        variableNamesOfDataChannels = [];
        independantVarRange = [];
        independantVarChannel = [];
        unitConversions = [];
        indiciesOfInterest = [];
        rawDataStructure
        remberSubsitutions 
        substitutions
        wheelPositionSuffix = [];
        returnVariablesWithWheelPositionSuffix 
        applyFilter = []; %[b,a] transfer function coeffs to send to filt filt if we are to filter data or [] for no filtering
        resampleDataToNPoints = [];
        returnOriginalMeas
        suppressWarnings       
    end
    
    methods
        function obj = jatecDaqDataSelector(varargin)
            
            %% Ohter defaults and substitutions
            %Deal with defaults
            persistent substitutions remberSubsitutionsFlag
            defaults = {'remberSubsitutions',[]};
            setDefaultsForVarargin(defaults,varargin);
            if isempty(remberSubsitutions)                                 %#ok did not mean to ref obj prop
                if ~isempty(remberSubsitutionsFlag)
                    remberSubsitutions = remberSubsitutionsFlag;           %#ok did not mean to ref obj prop
                end
            end
                
            if remberSubsitutions == false                                 %#ok did not mean to ref obj prop
                substitutions = [];
            else
                
            end
            obj.remberSubsitutions = remberSubsitutions;                   %#ok did not mean to ref obj prop
            obj.substitutions = substitutions;


            

        end %Creator method
        
        function obj = getChannelDataProcess(obj)
            %Very first thing, make sure that we actually have a rawData
            %field in teh daq file
            if ~isfield(obj.daq,'rawData')
                for iCh = 1:length(obj.channelsToGetData)
                    obj.rawDataStructure.(obj.variableNamesOfDataChannels{iCh}).meas = [];
                end
                warning('No rawData structure found on file: %s',displayDaqFiles(obj.daq,'suppressOutput',true))
                return
            end
            
            obj = setIndependantVariableRangesAndIndiciesOfInterest(obj);
            obj = seeIfChannelIsAvailableAndAddWheelPosIfNecessary(obj);
            if isempty(obj.channelsToGetData)
                fprinf('No channels to return, exiting\n')
                return
            end
            
            obj = getDataFromDaq(obj);
            obj = doubleCheckThatDataIsVectorsAndSameLength(obj);
            
            if ~isequal(obj.unitConversions,ones(size(obj.variableNamesOfDataChannels)))
                obj = applyUnitConversions(obj);
            end
            
            if ~isempty(obj.applyFilter)
                obj = filterData(obj);
            end
            
            if ~isempty(obj.resampleDataToNPoints)
                obj = resampleData(obj);
            end
            

        end%get the channel data
        
        function obj = setIndependantVariableRangesAndIndiciesOfInterest(obj)            
            %First special case where the appropiate indicies are passed in
            if strcmp(obj.independantVarChannel,'index')
                obj.indiciesOfInterest = obj.independantVarRange;
                return
            end
            
            %If it's empty grab all the data
            if isempty(obj.independantVarRange)
                %Grab any channel data
                channels = daqChannels(obj.daq);
                obj.indiciesOfInterest = 1:length(obj.daq.rawData.(channels{1}).meas);
                return;
            end                                
            
            %If not, we need to find the indexs in the right spot of the
            %data. Grab channel data and the range of it
            indepVariableChannelData = obj.daq.rawData.(obj.independantVarChannel).meas; 
            indepVariableRange = [min(indepVariableChannelData) max(indepVariableChannelData)];
                                        
            
            %See if it's a special case of a cell where we have some
            %special word
            if isa(obj.independantVarRange,'cell')
                for iEle = 1:length(obj.independantVarRange)
                    if strcmp(obj.independantVarRange{iEle},'end')
                        obj.independantVarRange{iEle} = indepVariableRange(2);
                    end
                end
                obj.independantVarRange = cell2mat(obj.independantVarRange);
            end
            
            %Now, it's gotta be a number array
            %A single data point
            if numel(obj.independantVarRange) == 1
                if obj.independantVarRange <= indepVariableRange(1) || obj.independantVarRange >= indepVariableRange(2) 
                    nearestPointInd = findNearestPoint(indepVariableRange,obj.independantVarRange);
                    if ~obj.suppressWarnings; warning('Requested data point is outside of available range, using nearest neigbor: %f',indepVariableRange(nearestPointInd));  end
                    obj.independantVarRange = indepVariableRange(nearestPointInd);
                end

            %A range of data, make sure start and finish are in the
            %range
            elseif numel(obj.independantVarRange) == 2             %range of points
                if obj.independantVarRange(1) < indepVariableRange(1)
                    if ~obj.suppressWarnings; warning('Requested data range outside of data, using starting value of %f',indepVariableRange(1)); end
                    obj.independantVarRange(1) = indepVariableRange(1);
                end

                if obj.independantVarRange(2) > indepVariableRange(2)
                    if ~obj.suppressWarnings; warning('Requested data range outside of data, using ending value of %f',indepVariableRange(2)); end
                    obj.independantVarRange(2) = indepVariableRange(2);
                end
                
            %Should must be a vector bigger than 1 turn that into a range
            else                                                   
                obj.independantVarRange = [min(obj.independantVarRange) max(obj.independantVarRange)];
            end

            %Now, get the indicies of interest
            if numel(obj.independantVarRange) == 1
                obj.indiciesOfInterest = findNearestPoint(indepVariableChannelData,obj.independantVarRange);
            else
                staringInd = findNearestPoint(indepVariableChannelData,obj.independantVarRange(1));
                endingInd  = findNearestPoint(indepVariableChannelData,obj.independantVarRange(2));
                obj.indiciesOfInterest = staringInd:endingInd;
            end
        end%setIndependantVariableRanges
        
        
        function obj = seeIfChannelIsAvailableAndAddWheelPosIfNecessary(obj)
        %The goal of this method is to see if a channel is available, and
        %if the user requested wheel positions, then augmnent the channel
        %list.
        
            %Deal with the wheel positions if passed in
            if ~isempty(obj.wheelPositionSuffix)
                if strcmp(obj.wheelPositionSuffix,'standard')
                    obj.wheelPositionSuffix = {'_L1','_R1','_L2','_R2'};
                end
            end
            
            allChannels = daqChannels(obj.daq);
            deleteChannel = [];
            %If these are wheel positions
            for iCh = 1:length(obj.channelsToGetData)
                channel = obj.channelsToGetData{iCh};
                variableName = obj.variableNamesOfDataChannels{iCh};
                conversion = obj.unitConversions(iCh);
                
                if ~isempty(obj.wheelPositionSuffix)
                    for iPos = 1:length(obj.wheelPositionSuffix)
                        tempChannel = [channel obj.wheelPositionSuffix{iPos}];
                        if obj.returnVariablesWithWheelPositionSuffix
                            tempVariableNameOfData = [variableName obj.wheelPositionSuffix{iPos}];
                        else
                            tempVariableNameOfData = variableName;
                        end
                            
                        

                        if ismember(tempChannel,allChannels)
                            %See if we're appending the wheel position
                            %suffix
                            obj.variableNamesOfDataChannels = [obj.variableNamesOfDataChannels; tempVariableNameOfData];
                            obj.variableNamesOfDataChannels = obj.variableNamesOfDataChannels;
                            obj.channelsToGetData = [obj.channelsToGetData; tempChannel];
                            obj.rawDataStructure.(tempVariableNameOfData).channelName = tempChannel;
                            obj.unitConversions = [obj.unitConversions; conversion];
                        else
                            if ~obj.suppressWarnings; warning('Channel %s not found!',tempChannel); end
                        end
                    end
                    deleteChannel = [deleteChannel; {channel}];
                else
                    
                    %If not wheel positions, see if the channel exist
                    if ~ismember(channel,[allChannels; 'index'])
                        obj.remberSubsitutions = false;
                        askUserForNewChannel = false;
                        newChannelSelected = false;

                        %This is the first time
                        if isempty(obj.remberSubsitutions)
                            choice = questdlg('Would you like to remember subsitutions?','Track subsitutions',...
                                              'Yes','No','Yes');
                            switch choice
                                case 'Yes'
                                    obj.remberSubsitutions = true;
                                case 'No'
                                    obj.remberSubsitutions = false;
                            end
                        end

                        %See if subsitutions exist
                        if ~isempty(obj.substitutions)
                            [subExists,ind] = ismember(channel,obj.substitutions(:,1));
                            if subExists
                                newChannelName = obj.substitutions{ind,2};
                                if ismember(newChannelName,daqChannels(obj.daq));
                                    askUserForNewChannel = false;
                                end
                            end
                        end

                        if askUserForNewChannel
                            [newChannelName,newChannelSelected] = daqChannelFilter(obj.daq,'figTitle',sprintf('Choose replacement channel for: %s',channel),'addCancelButton',true,'cancelButtonText','Skip Channel','multiSelect',false);
                            if obj.remberSubsitutions
                                newSub = {channel newChannelName};
                                obj.substitutions = [obj.substitutions; newSub];
                            end
                        end


                        if ~newChannelSelected
                            obj.rawDataStructure.(obj.variableNamesOfDataChannels{iCh}).channelName = [];
                            obj.daq.rawData.(channel).meas = [];
                        else
                            obj.rawDataStructure.(obj.variableNamesOfDataChannels{iCh}).channelName = newChannelName;

                            %Ask user to duplicate the channel name or rename
                            %the channel on the daq file
                            tempChoice = sprintf('Neither, temporarily use "%s"',newChannelName);
                            dlgTxt =  {sprintf('Rename channel: "%s" as "%s"',newChannelName,channel)
                                       sprintf('Or, duplicate "%s" as "%s" for next time',newChannelName,channel) };
                            choice = questdlg(dlgTxt, 'Rename/duplicate channel', ...
                                'Rename','Duplicate',tempChoice,'Rename');
                            % Handle response
                            switch choice
                                case 'Rename'
                                    obj.daq = renameDaqChannel(newChannelName,channel,obj.daq);
                                case 'Duplicate'
                                    obj.daq = duplicateChannel(newChannelName,channel,obj.daq);
                            end

                        end
                    else
                        obj.rawDataStructure.(variableName).channelName = channel;
                    end
                end
            end
            
            %Delete anychannels necessary
            ind = [];
            for iCh = 1:length(deleteChannel)
                [~,trueFalse] = ismember(obj.channelsToGetData,deleteChannel{iCh});
                ind = [ind find(trueFalse)];
            end
            obj.channelsToGetData(ind) = [];
            obj.variableNamesOfDataChannels(ind) = [];
            obj.unitConversions(ind) = [];
            
            %Last thing check that there is a daq structure available for
            %the data channel
            deleteChannel = {};
            for iCh = 1:length(obj.channelsToGetData)
                channel = obj.channelsToGetData{iCh};
                if ~strcmp(channel,'index')
                    if ~isa(obj.daq.rawData.(channel),'struct')
                        if ~obj.suppressWarnings; warning('Channel %s is not a daq strucutre and cannot be returned',channel); end
                        deleteChannel = [deleteChannel; channel];
                    end
                end
            end
            
            %Delete anychannels necessary
            ind = [];
            for iCh = 1:length(deleteChannel)
                obj.rawDataStructure = rmfield(obj.rawDataStructure,deleteChannel{iCh});
                [~,trueFalse] = ismember(obj.channelsToGetData,deleteChannel{iCh});
                ind = [ind find(trueFalse)];
            end
            obj.channelsToGetData(ind) = [];
            obj.variableNamesOfDataChannels(ind) = [];
            obj.unitConversions(ind) = [];            
            

        end%seeIfChannelIsAvailable
        
        function obj = getDataFromDaq(obj)
            for iCh = 1:length(obj.channelsToGetData)
                channel = obj.rawDataStructure.(obj.variableNamesOfDataChannels{iCh}).channelName;   
                if strcmp(channel,'index')
                    daqChannelStructure.name = 'index';
                    daqChannelStructure.units = '';
                    daqChannelStructure.mathChannel = true;
                    meas = [];
                    %Try to get a monoticially increasing channel
                    availChannels = daqChannels(obj.daq);
                    for iChAvail = 1:length(availChannels)
                        if ismonotonic(obj.daq.rawData.(availChannels{iChAvail}).meas)
                            meas = 1:length(obj.daq.rawData.(availChannels{iChAvail}).meas);
                            break
                        end
                    end
                    %Just in case there wasn't a monotically increasing
                    %channel, just grab the first one
                    if isempty(meas)
                        meas = 1:length(obj.daq.rawData.(availChannels{1}).meas);
                    end
                    daqChannelStructure.meas = rowVector(meas);
                    
                    
                    obj.rawDataStructure.(obj.variableNamesOfDataChannels{iCh})             = catstruct(obj.rawDataStructure.(obj.variableNamesOfDataChannels{iCh}),daqChannelStructure);
                elseif ~isempty(channel)
                    daqChannelStructure = obj.daq.rawData.(channel);
                    obj.rawDataStructure.(obj.variableNamesOfDataChannels{iCh})             = catstruct(obj.rawDataStructure.(obj.variableNamesOfDataChannels{iCh}),daqChannelStructure);
                   
                    if ~isfield(obj.rawDataStructure.(obj.variableNamesOfDataChannels{iCh}),'meas')
                        obj.rawDataStructure.(obj.variableNamesOfDataChannels{iCh}).meas        = [];
                    else
                        %first make sure it's going to be a double so we
                        %can get the indicies, if not, then just return the
                        %meas field what ever it is
                        
                        %If the user wants the originalMeas, see if it is there
                        if obj.returnOriginalMeas && isfield(obj.rawDataStructure.(obj.variableNamesOfDataChannels{iCh}),'originalMeas')
                            channelMeasField = obj.rawDataStructure.(obj.variableNamesOfDataChannels{iCh}).originalMeas;
                        else
                            channelMeasField = obj.rawDataStructure.(obj.variableNamesOfDataChannels{iCh}).meas;
                        end
                        
                        
                        
                        if isa(channelMeasField,'double') && isempty(find(~ismember(obj.indiciesOfInterest,1:length(channelMeasField)), 1))
                            obj.rawDataStructure.(obj.variableNamesOfDataChannels{iCh}).meas = channelMeasField(obj.indiciesOfInterest);
                        else
                            obj.rawDataStructure.(obj.variableNamesOfDataChannels{iCh}).meas = channelMeasField;
                        end
                        
                        
                    end
                else
                    obj.rawDataStructure.(obj.variableNamesOfDataChannels{iCh}).meas = [];
                end
                
                if ~isfield(obj.rawDataStructure.(obj.variableNamesOfDataChannels{iCh}),'name')
                    obj.rawDataStructure.(obj.variableNamesOfDataChannels{iCh}).name = channel;
                end
                
                if ~isfield(obj.rawDataStructure.(obj.variableNamesOfDataChannels{iCh}),'units')
                    obj.rawDataStructure.(obj.variableNamesOfDataChannels{iCh}).units = channel;
                end
                
                    
            end
        end%getDataFromDaq        
        
        function obj = doubleCheckThatDataIsVectorsAndSameLength(obj)
            channels = fieldnames(obj.rawDataStructure);
            for iCh = 1:length(channels)
                channel = channels{iCh};
                daqChannelName{iCh} = obj.rawDataStructure.(channel).channelName;
                meas = obj.rawDataStructure.(channel).meas;
                [r,c] = size(meas);
                if r == 1 || c == 1
                    if r ==1 %If not it's a matrix and don't trasnpose
                        meas = rowVector(meas);
                        obj.rawDataStructure.(channel).meas = meas;
                        obj.daq.rawData.(daqChannelName{iCh}).meas = meas;
                        if ~obj.suppressWarnings; warning('Channel %s on daq file transposed for correct dimensions',daqChannelName{iCh}); end
                        
                    end
                end
                
                [r,c] = size(meas);
                nRows(iCh) = r;
                nCols(iCh) = c;
            end
                
            [~,ind] = unique(nRows);
            if numel(ind) > 1
                if ~obj.suppressWarnings; 
                    displayDaqDataSizesAndStats(obj.daq);
                    warning('Number of measurments not consisent in the data file'); 
                end
            end
        end %obj = doubleCheckThatDataIsVectorsAndSameLength(obj)
        
        function obj = filterData(obj)
            %This funciton will apply a passed in coeffs to filtfilt
            b = obj.applyFilter{1};
            a = obj.applyFilter{2};
            for iCh = 1:length(obj.variableNamesOfDataChannels)
                daqChannel = obj.channelsToGetData{iCh};
                channel = obj.variableNamesOfDataChannels{iCh};
                
                filtData = filtfilt(b,a,obj.rawDataStructure.(channel).meas);
                obj.rawDataStructure.(channel).meas = filtData;
                
                %see if original data exist, if not save it before
                %modifying
                if ~isfield(obj.daq.rawData.(daqChannel),'originalMeas')
                    obj.daq.rawData.(daqChannel).originalMeas = obj.daq.rawData.(daqChannel).meas;
                end
                obj.daq.rawData.(daqChannel).meas = filtData;
            end
        end%Filter data
        
        function obj = resampleData(obj)
            %This funciton will resample the data to number of requested
            %points
        end%resample data
        
        function obj = applyUnitConversions(obj)
            for iCh = 1:length(obj.variableNamesOfDataChannels)
                daqChannel           = obj.channelsToGetData{iCh};
                channelToBeAssigned  = obj.variableNamesOfDataChannels{iCh};
                
                oldData = obj.rawDataStructure.(channelToBeAssigned).meas;
                newData = oldData*obj.unitConversions(iCh);

                obj.daq.rawData.(daqChannel).meas = newData;
                obj.rawDataStructure.(channelToBeAssigned).meas = newData;
            end
        end
    end
end