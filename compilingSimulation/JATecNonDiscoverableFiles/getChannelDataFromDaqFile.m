function [daq, rawDataStruct] = getChannelDataFromDaqFile(daq,vars,varargin)
%This function returns variables for all daq channels that you define in
%the vars cell array.  Vars is dimenstion nx2 or nx3 where (:,1) is the
%name of the variable that you would like it defined as and (:,2) is the
%actual name as it appears in the daq structure. The third colum is used to
%convert units of the returned variables (this is optional). For example:
%   vars = {'t' 'time'
%           'vx' 'Vx'};
%
%   vars = {'t'   'time'  1
%            'vx' 'Vx'    myConstants.kph2mps};
%will return a variables t and vx in the workspace of the funciton calling
%this function.  It will call them like: t = daq.rawData.(vars{i,2}).meas
%*vars{i,3}. Therefore daq is the local instance of a single daq file from
%the global DAQ cell array
%
%This was written to help deal with all the variaion in channel naming.
%
%INPUTS:
%  daq  -  instance of a daq file                              class struct
%  vars -  (optional) cell array {'var name you want' 'daq channel name'}   
%          or just a row array of channels
%          or 'gui' to select                            class cell or char
%  varargin - used to set default values                                              
%          
% 
%OUTPUTS:
%  daq - updated daq file NOT just the requested time          class struct
%  daqDataStruct - this is mainly used for making axes labels.  This is a
%        daq.rawData style output in the following format:
%        rawDataSruct.(vars{iRow,1}) = daq.rawData.(vars{iRow,2}).
%        This is updated to be the requested time range.
%
%Created: 20 Jan 2013 - Jeff Anderson 
%Updated: 19 Jun 2014 - Jeff Anderson - included persistant variable for
%                       channels skipped so when looping through several
%                       files.  As well as a persistent var for keeping up
%                       with substitutions
%Updated: 03 Mar 2015 - Jeff Anderson, flipped input so that you could
%                       leave off vars and it would default to all
%                       channels.  Also made it so you still use legacy
%                       code and it'd fix the inputs
%Updated: 27 Mar 2015 - JRA - have it return [] vars if you click skip, most
%                       function calling this will expect a var with this 
%                       name
%Updated: 02 Oct 2015 - JRA - error out if user closes selection dialog
%Updated: 04 Oct 2015 - JRA - save new channel selections 
%Updated: 04 Oct 2015 - JRA - new output "rawDataStruct" to return a daq
%                       style structure of the vars requested for ease of 
%                       making labels on axes.
%Updated: 14 Dec 2015 - JRA - changed "rawDataStruct" to have a name if
%                       none is there, used for plotting
%Updated: 11 Feb 2016 - jra - added a 'Name' to the rawDataStructure, for
%                       the case where the channel only has a meas field,
%                       it will then at least return that channel name.
%Updated: 12 May 2016 - jra - made sure units came back in the raw data
%                       structure. Used ' ' if no units in daq
%Updated: 04 Nov 2016 - jra - Made it so you can pass in a char
%Updated: 05 Jan 2017 - jra - dramatic departure and using
%                       jatecDaqDataSelector object to select the data
%Updated: 16 Jan 2017 - jra - moved options to varargin
%Updated: 13 Feb 2017 - jra - option to remove nans
%Updated: 14 May 2017 - jra - filtering option, unit conversion added to
%   the vars input. 'wheelPositionSuffix' option to return variables w/o
%   wheel position suffix
%Updated: 11 Sep 2017 - jra - added 'index' as a special channel to get.
%Updated: 10 Oct 2017 - Jeff Anderson - added option to return origMeas
%    from the daq file.
%Updated: 10 Oct 2017 - Jeff Anderson - origMeas changed to originalMeas 
%    per documentation
%Updated: 23 Feb 2018 - Jeff Anderson. Bug fix, for the case of an invalid
%    file, it won't error out; gracefully exits now.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% Varargin
defaults = {'indepVariableChannel',[]  %What channel to look for data in can be directly 'index' to specify data indicies. Use 'atIndepVariableTime' to specify indicies
            'atIndepVariableTime', []  %What range of data to look in
            'wheelPositionSuffix',[]
            'returnVariablesWithWheelPositionSuffix',true % If wheelPosition suffix isn't empty it will return variables for each channel, if this option if false it will return fy insteady of fy_L1 etc.. CAREFUL! It will override for multiple of the same channels
            'removeNans',false
            'applyFilter',[]
            'showWarningForNanRemoved', false
            'resampleDataToNPoints',false,
            'returnOriginalMeas',false
            'suppressWarnings',false
            'returnUpdatedDaq',false};  %If we filter or unit convert, the file changes, return the udpated files??
setDefaultsForVarargin(defaults,varargin)

%Option Definitions
%  atIndepVariableTime - used to select a specific point in time/dist or
%          range of time/dist. Pass in either a single double value for a
%          specific point or a [startTime stopTime] vectory for a range of
%          points. If no argument is used, then the full time history is
%          returned.                                    class double arrary
%  indepVariableChannel - if 'atIndepVariableTime' is used, then what
%          channel is the independant varaible? Usually time or distance,
%          if no argument is used, then the function
%          getIndependantVariableChannel() will try to find the correct
%          one.    
%  wheelPositionSuffix - pass in cell array of wheel positions to get all
%          the data from those channels.  Alternately pass in 'standard' to
%          use the standard {'L1','R1','L2','R2'} wheel positions.
% removeNans - This will remove all nans
% applyFilter - pass in {b,a} and it will apply filtfilt to each channel.
%          Must be a cell to pass in multiple vectors of coeffs.
% resampleDataToNPoints = will resample the data down to the number of
%          points requested. Useful for cross plots.


%% Deal w/ input args

%I flipped the inputs so it would be daq, then vars not vice versa.  This
%matches convention used in toher files.  So need to take care of it if it
%is backwards from legacy code.
if nargin == 2
    if isa(vars,'struct') && isa(daq,'cell')
        tempDaq = vars;
        tempVars = daq;
        daq = tempDaq;
        vars = tempVars;
    end
end

%Need to see if the first loaded daq is the intended file, if so grab it
loadFirstDaq = false;
if ~exist('daq','var'); 
    loadFirstDaq = true;
elseif isempty(daq)
    loadFirstDaq = true;
end

useGlobal = false;
if loadFirstDaq
    global DAQ;  %#ok, don't need to global unless we're loading this file
    daq = DAQ{1}; 
    useGlobal = true;
end


%Deal with vars argument, if there is none, then assume all channels
if ~exist('vars','var'); 
    vars = [daqChannels(daq) daqChannels(daq)]; 
elseif isempty(vars)
    vars = [daqChannels(daq) daqChannels(daq)]; 
end
    



%If it's a char, see if it's 'gui', if so prompt user, if not, then it's
%just one channel
if isa(vars,'char')
    %First see if it's 'gui' then pompt user
    if strcmp(vars,'gui')
        channels = daqChannelFilter(daq);
        vars = [channels channels];
    else
        temp = vars;
        clear vars;
        vars{1} = temp;
    end
end


%If a cell array vector is passed in, assume the channel names are the
%variables to be assigned.
% [r,c] = size(vars);
% if isvector(vars) && c == 1
%     vars = rowVector(vars);
%     vars = [vars vars]; 
% elseif numel(vars)
    [~,c] = size(vars);
    if c == 1
        vars = [vars vars]; 
    end
% end



%see if there are any duplicate channels requested. Make sure we're
%assiging a unique set of variables.
for iCol = 1:1 
    [~,ind] = unique(vars(:,iCol));
    if length(ind) ~= length(vars(:,iCol))
        error('Requested channel names must be a unique set!')
    end
end

%Conversions
[r,c] = size(vars);
if c == 2 %Assume no conversions
    for iCh = 1:r
        vars{iCh,3} = 1;
    end
end

%Wheel position suffix: make sure it's a cell array
if isa(wheelPositionSuffix,'char')
    temp = wheelPositionSuffix;
    clear wheelPositionSuffix
    wheelPositionSuffix{1} = temp;
end


%% Meat and potatos



%Get the data
dataSelector = jatecDaqDataSelector;
dataSelector.daq = daq;
dataSelector.channelsToGetData = vars(:,2);
dataSelector.unitConversions   = cell2mat(vars(:,3));
dataSelector.returnOriginalMeas = returnOriginalMeas;
dataSelector.variableNamesOfDataChannels = vars(:,1);
dataSelector.independantVarRange = atIndepVariableTime;
dataSelector.independantVarChannel = indepVariableChannel;
dataSelector.wheelPositionSuffix = wheelPositionSuffix;
dataSelector.returnVariablesWithWheelPositionSuffix = returnVariablesWithWheelPositionSuffix;
dataSelector.applyFilter = applyFilter;
dataSelector.suppressWarnings = suppressWarnings;

dataSelector.resampleDataToNPoints = [];

dataSelector = dataSelector.getChannelDataProcess();

if returnUpdatedDaq
    daq = dataSelector.daq; %if we filter or change the units, the file changes
end


% obj = jatecDaqDataSelector(daq,vars,atIndepVariableTime,indepVariableChannel);
% daq = obj.daq;
rawDataStruct = dataSelector.rawDataStructure;
if isempty(rawDataStruct)
    return
end

channels = fieldnames(rawDataStruct);

%Remove nans?
ind = [];
if removeNans
    for iCh = 1:length(channels)
        channel = channels{iCh};
        meas = rawDataStruct.(channel).meas;
        ind = [ind; find(isnan(meas))];
    end
    ind = unique(ind);
    if ~isempty(ind)
        if showWarningForNanRemoved
            warning('Nan indicies removed!')
            disp(ind)
        end
    end
end


for iCh = 1:length(channels)
    channel = channels{iCh};
    meas = rawDataStruct.(channel).meas;
    meas(ind) = [];
    assignin('caller',channel,meas)
end

if useGlobal
    DAQ{1} = daq;
end



%%% OLD DEPRICATED CODE %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% persistent skipChannel
% persistent subChannel %Var to store substitutions subChannel{:,1} = old name, subChannel{:,2} = new name
% 
% 
% 
% %% Get them from daq files
% [nCh,nCols] = size(vars);
% if nCols == 2 && nCh == 1
%     warning('should vars be a row vector? only channel %s returned as %s',vars{1,2},vars{1,1});
% end
% 
% if nCols == 1 %Need to make the assigned variabled and channel name the same:
%     vars = [vars vars];
% end
% for iCh = 1:nCh
%     daqChannel = vars{iCh,2};
%     if isfield(daq.rawData,daqChannel)                                     %Check for the channel
%         if isfield(daq.rawData.(daqChannel),'meas')                         %Some times measurements aren't a field
%             %Get the daq channel data and save this in the rawDataStruct
%             %output
%             temp = daq.rawData.(daqChannel).meas;
%             
%             
%             if isa(temp,'double')
%                 temp = rowVector(temp);
%                 assignin('caller',vars{iCh,1},temp);
%                 rawDataStruct.(vars{iCh,1}) = daq.rawData.(vars{iCh,2});   %Grab teh full structure
%                 
%                 %Check if there is a name or not?
%                 if ~isfield(rawDataStruct.(vars{iCh,1}),'name') 
%                     rawDataStruct.(vars{iCh,1}).name = vars{iCh,2};
%                 elseif  strcmp(rawDataStruct.(vars{iCh,1}).name,' ')
%                     rawDataStruct.(vars{iCh,1}).name = vars{iCh,2};
%                 end
%                 
%                 %Check if there are units or not, if not return ' '
%                 if ~isfield(rawDataStruct.(vars{iCh,1}),'units') 
%                     rawDataStruct.(vars{iCh,1}).units = ' ';
%                 end
%                 
%             else
%                 warning('Channel: %s skipped. Identified as class %s and needs to be class double on file: %s',daqChannel,class(temp),daq.header.filename)
%             end
%         end
% 
%     else %channel is not found
%         proptToSelectNewChannel = true; %flag to see if we need to prompt user for a new channel selection
%         
%         %First see if we have already skipped this
%         try ind = find(ismember(skipChannel,daqChannel));
%             if ~isempty(ind)
%                proptToSelectNewChannel = false;
%                warning('Will skip %s again, this var was previously skipped',daqChannel);
%                temp = [];%nan(length(daq.rawData.time.meas),1);
%                assignin('caller',vars{iCh,1},temp)
%                rawDataStruct.(vars{iCh,1}) = createDaqChannelData([],' ',sprintf('Channel %s not found on file %s',vars{iCh,2},daq.header.filename));
%             end
%         end
%         
%         %If not, then see if we already hav eit assigned with a substitution
%         if ~isempty(subChannel)
%             try ind = find(ismember(subChannel(:,1),daqChannel));
%                 if ~isempty(ind)
%                     proptToSelectNewChannel = false;
%                     warning('Variable %s was assigned with daq channel %s instead of %s',vars{iCh,1},subChannel{ind,2},daqChannel);
%                     ind = find(strcmp(subChannel{:,1},daqChannel));
%                     temp = daq.rawData.(subChannel{ind,2}).meas;
%                     assignin('caller',vars{iCh,1},temp);
%                     rawDataStruct.(vars{iCh,1}) = daq.rawData.(subChannel{ind,2});
%                     if ~isfield(rawDataStruct.(vars{iCh,1}),'name')
%                         rawDataStruct.(vars{iCh,1}).name = vars{iCh,2};
%                     end
%                 end
%             end
%         end        
%         
%             
%             
%         if proptToSelectNewChannel    
%             allChannels = fieldnames(daq.rawData);
%             prompString = sprintf('Channel: %s is not found, choose replacement:',daqChannel);
%             warning(prompString);                                          %I did this incase operating in command mode, not sure what channel errored out
%             %Make a figure for selection (can't use listdlg bec we want
%             %skip and cancel
%             hPrompFigure = figure('Name','Channel Replacement','Color','w','ToolBar','none','MenuBar','none');
%             vSpace = uiextras.VBox('Parent',gcf,'Padding',10);
%             
%             %What channel
%             uiextras.Empty('Parent',vSpace,'Visible','off');
%             uicontrol('Parent',vSpace,'Style','Text','String',prompString,'BackgroundColor','w');
%             hList = uicontrol('Parent',vSpace,'Style','list','string',allChannels,'BackgroundColor','w');
%             
% 
%             
%             %Buttons
%             hButtonSpace = uiextras.HBox('Parent',vSpace);
%             uicontrol('Parent',hButtonSpace,'String','OK',          'callback','set(gcf,''UserData'',''ok'',''visible'',''off'')')
%             uicontrol('Parent',hButtonSpace,'String','Skip Channel','callback','set(gcf,''UserData'',''skip'',''visible'',''off'')')
%             uicontrol('Parent',hButtonSpace,'String','Cancel',      'callback','set(gcf,''UserData'',''cancel'',''visible'',''off'')')
%             
%             %Fix uiextras spaces
%             set(vSpace,'Sizes',[20 40 -1 50])
%             
%             %Callback for being delted           
%             waitfor(gcf,'visible','off')
%             selection = get(gcf,'userdata');
%             iSelection = get(hList,'Value');
%             delete(gcf)
%             
%             switch selection
%                 case 'ok'
%                     %Assign the var
%                     newChannelName = allChannels{iSelection};
%                     temp = daq.rawData.(newChannelName).meas;
%                     assignin('caller',vars{iCh,1},temp);
%                     rawDataStruct.(vars{iCh,1}) = daq.rawData.(newChannelName);
%                     if ~isfield(rawDataStruct.(vars{iCh,1}),'name')
%                         rawDataStruct.(vars{iCh,1}).name = vars{iCh,2};
%                     end
%                     
%                     %Ask the user if they want to change the var name
%                     questString = sprintf('Would you like to replace the channel on the daq file from %s to %s?',allChannels{iSelection},vars{iCh,2});
%                     choice = questdlg(questString, ...
%                                       'Rename daq channel', ...
%                                       'Yes','No','Yes');
%                     if strcmp(choice,'Yes')                        
%                         daq.rawData.(vars{iCh,2}) = daq.rawData.(allChannels{iSelection});
%                         daq.rawData = rmfield(daq.rawData,allChannels{iSelection});
%                     end
% 
%                     %Note the substitution
%                     subChannel{length(subChannel)+1,1} = daqChannel;
%                     subChannel{length(subChannel),2} = allChannels{iSelection};
%                 case 'skip'
%                     skipChannel{length(skipChannel)+1} = daqChannel;
%                     temp = []; % Fix jra 8/12/15 % Return [] instead of nan nan(length(daq.rawData.time.meas),1);
%                     assignin ('caller',vars{iCh,1},temp)
%                     rawDataStruct.(vars{iCh,1}) = createDaqChannelData([],' ',sprintf('Channel %s not found on file %s',vars{iCh,2},daq.header.filename));
%                 case 'cancel'
%                     error('Process aborted')
%             end
%         end %Promp user to select a new channel
%     end %if it is a channel or not
%     
% end %loop through channels
% 
% % %The last thing is to return a special time or set of times
% % if ~isempty(indepVariableChannel) && ~isempty(atIndepVariableTime)
% %     disp('sac')
% %     channels = fieldnames(rawDataStruct)
% %     for iCh = 1:length(fields
% %     rawDataStruct.
% %     
% % end
% 
% 
% end%getchanneldatafromdaqfile
% 
% 
