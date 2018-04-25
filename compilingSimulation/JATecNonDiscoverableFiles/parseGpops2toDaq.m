function daq = parseGpops2toDaq(output,states,controls,indepVar,units,names,varargin)
%This function will parse the GPOPS-II optimal control software code to a
%daq file.
%INPUT:
%    output - GPOPS-II output structure                        class struct
%    states - cell array of channel names for the states         class cell
%    control - cell array of channel names for controls          class cell
%    indepVar - char of the channel name of the independance var class char
%    units - (optional)  cell array of units for states, control, indep in 
%            that order                                          class cell
%    varargin - optional - used to set defaults below.
%
%Creation: 5 May 2015 - Jeff Anderson
%Update: 29 Sep 2015 - Jeff Anderson - removed teh frenet calcs here
%Update: 04 Mar 2016 - Jeff Anderson - added names
%Update: 06 Feb 2017 - Jeff Anderson
%    If channel is standard, use addMathChannelsThatAreStandardChannels
%    rather than createDaqChannelData to add data.
%Updated: 27 Nov 2017 - Jeff Anderson - grab a different field of the
%    solution, the interpsolution.
%Updated: 05 Feb 2018 Jeff Anderson - minor fixes so it didn't
%    automatically look for standard channel names.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

defaults = {'solutionField' 'interpsolution' %field in the gpops results output.results.(solutionField) originally just solution but looks like the interpsol is more what I need.
            'lookForStandardChannelNames' false
            'reinterpolateDataToIndependantVariable',true
            'dIndependantVariable',0.25};
setDefaultsForVarargin(defaults,varargin);


%% Input error checking
if nargin <4
    error('Expecing at least 4 inputs, only %i called',nargin)
elseif nargin<5
    units = cell(length(states)+length(controls)+1,1);
end


%Check states input
[~,nStates] = size(output.result.(solutionField).phase(1).state);
[rState,c] = size(states);
if ~isa(states,'cell')
    error('State vector not a cell array')
elseif c > 1
    error('State vector should be row vector')
elseif nStates ~= rState
    error('Incorrect number of states defined, expected %i, received %i',nStates,rState)
end

%Check control

if isfield(output.result.(solutionField).phase(1),'control')
    [~,nControl] = size(output.result.(solutionField).phase(1).control);
    [rControl,c] = size(controls);
    if ~isa(controls,'cell')
        error('Control vector not a cell array')
    elseif c > 1
        error('Control vector should be row vector')
    elseif nStates ~= rState
        error('Incorrect number of controls defined, expected %i, received %i',nControl,rControl)
    end
else
    nControl = 0;
end

%Check indepVar
if ~isa(indepVar,'char')
    error('Expected char, received %s',class(indepVar))
end

%Check units
[rUnits,c] = size(units);
if c > 1
     error('Units vector should be row vector')
elseif rUnits ~= nStates+nControl+1
    error('Units should be row cell vector with %i rows, recieved %i rows',nStates+nControl+1,rUnits)
end

if ~exist('names','var')
    names = cell(length(states) + nControl + 1);
end

%% Put the channels names in a struct so we can get to them later
daq.gpopsNames.states   = states;
daq.gpopsNames.controls = controls;
daq.gpopsNames.indepVar = indepVar;
daq.gpopsNames.units = units;
daq.gpopsNames.names = names;

%% Parsing code



daq.gpopsOutput = output;
name = '';
try name = output.name; end
    
daq.header.filename = sprintf('%s GPOPS Output - %s',datestr(now,'yyyy-mm-dd HH:MM'),name);
daq.header.path = pwd;
daq.header.creationDate.datenum = datenum(now);
daq.header.creationDate.day = datestr(now,'DD');
daq.header.creationDate.month = datestr(now,'MM');
daq.header.creationDate.year = datestr(now,'YYYY');
daq.header.creationDate.hms = datestr(now,'HH:MM:ss');

if ~isfield(daq,'rawData')
    daq.rawData = [];
end


%Grab teh actual data
for iPhase = 1:length(output.result.(solutionField).phase)
    if iPhase == 1
        channelCount = 1;
        for iState = 1:length(states)
            stateData = output.result.(solutionField).phase(iPhase).state(:,iState);
            stateName = states{iState};
            if isStandardChannel(stateName) && ~lookForStandardChannelNames
                daq.rawData = addMathChannelsThatAreStandardChannels(daq.rawData,stateName,stateData);
            else
                daq.rawData.(stateName) = createDaqChannelData(stateData,units{channelCount},names{channelCount});
            end
            channelCount = channelCount+1;
        end

        for iControl = 1:nControl
            controlName = controls{iControl}; 
            controlData = output.result.(solutionField).phase(iPhase).control(:,iControl);
            if isStandardChannel(controlName) && ~lookForStandardChannelNames
                daq.rawData = addMathChannelsThatAreStandardChannels(daq.rawData,controlName,controlData);
            else
                daq.rawData.(controlName) = createDaqChannelData(controlData,units{channelCount},names{channelCount});
            end
            channelCount = channelCount+1;
        end

        indepSol = output.result.(solutionField).phase(iPhase).time;
        if isStandardChannel(indepVar) && ~lookForStandardChannelNames
            daq.rawData = addMathChannelsThatAreStandardChannels(daq.rawData,indepVar,indepSol);
        else
            daq.rawData.(indepVar) = createDaqChannelData(indepSol,units{channelCount},names{channelCount});
        end
        
    else
        for iState = 1:length(states)
            solState = output.result.(solutionField).phase(iPhase).state(:,iState);
            stateData = [daq.rawData.(states{iState}).meas; solState];
            stateName = states{iState};
            daq.rawData.(stateName).meas = stateData;
        end

        for iControl = 1:nControl
            controlName = controls{iControl};
            controlData = [daq.rawData.(controls{iControl}).meas; contSol];
            contSol = output.result.(solutionField).phase(iPhase).control(:,iControl);
            daq.rawData.(controlName).meas = controlData;
        end

        indepSol = output.result.(solutionField).phase(iPhase).time;
        daq.rawData.(indepVar).meas = [daq.rawData.(indepVar).meas; indepSol];
        
    end
    

end

%% Re-interpolate because not gauranteed to be on a unifrom grid from gpops
if reinterpolateDataToIndependantVariable
    oldIndepVar = daq.rawData.(indepVar).meas;
    newIndepVar = rowVector(daq.rawData.(indepVar).meas(1):dIndependantVariable:daq.rawData.(indepVar).meas(end));
    channels = daqChannels(daq);
    for iCh = 1:length(channels)
        channel = channels{iCh};
        daq.rawData.(channel).meas = interp1(oldIndepVar,daq.rawData.(channel).meas,newIndepVar);
    end
end

% %% Frenet frame calcs
% if ~isempty(track)
% else
%     error('Expected track name to be loaded as global "track" variable');
% end
% 
% s = daq.rawData.distance.meas;
% trackData = trackParametersAtLoopedDistance(track,s);
% k = trackData.curvature.meas;
% trackWidth = trackData.trackWidth.meas(1);
% 
% rawData2 = frenetToXYFrameTrajectory(k,s,daq.rawData.ey.meas,...
%                                          daq.rawData.ePsi.meas,...
%                                      pathX0,pathY0,pathPsi0,...
%                                      trackWidth,'suppressPlot',true);
%                                  
% daq.rawData = catstruct(daq.rawData,rawData2);




