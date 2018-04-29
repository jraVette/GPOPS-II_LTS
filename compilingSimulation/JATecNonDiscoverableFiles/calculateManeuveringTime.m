function [times, daq] = calculateManeuveringTime(daq,varargin)
%This function will calculate the maneuvering time for the optimal control
%simulations
%
%INPUTS:
%    daq - optional local instance of the daq file             class struct
%    varargin - used to set 
%OUTPUT:
%    times - lap times of the daq files                        class double
%    daq - updated local instance of the file                  class struct
%
%Creation: 8 Mar 2016 - Jeff Anderson
%Updated:  3 Jan 2017 - Jeff Anderson: error checking to be sure that the
%    simulation acutally made it to the finish line.
%Updated:  9 Jan 2017 - Jeff Anderson: changed all track to have
%    startDistance and finishDistance fields.  This should eliminate
%    confusions.
%Updated: 18 Aug 2017 - Jeff Anderson - made the fields in
%    generateInitialDaq.m provide the start/stop if available.  Needed this
%    for short segments.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

if ~exist('daq','var'); daq = []; end
DAQ = useLocalOrGlobalDaqFiles(daq);

defaults = {'suppressOutput',false}; 
setDefaultsForVarargin(defaults,varargin)


for iFile = 1:length(DAQ)
    daq = DAQ{iFile};
    distance = daq.rawData.distance.meas;
    
    %Grab the track starting distance
    if isfield(daq.header,'timingDistanceStart')
        startDistance = daq.header.timingDistanceStart;                    %Set up in generateInitialdaq
    elseif isfield(daq.header.track,'startDistance')
        startDistance = daq.header.track.startDistance;
    else
        startDistance = 0;
        warning('No daq.header.timingDistanceStart found!!! assuming s0 = 0m');
    end
    
    %Get the finish distance
    if isfield(daq.header,'timingDistanceFinish')
        endDist = daq.header.timingDistanceFinish;
    elseif isfield(daq.header.track,'finishDistance')
        endDist = daq.header.track.finishDistance;
    else
        warning('No daq.header.timingDistanceFinish found!!! Using track end distance');
        endDist = daq.header.track.distance.meas(end);
    end
    
    
    %Get start/finish indicies
    indStart = findNearestPoint(distance,startDistance);
    indEnd = findNearestPoint(distance,endDist);
    lapTime = nan;
    
    %Do some error checking and make sure that the last distance recorded
    %is indeed > than the endDistnace
    if max(distance) < endDist
        warning('Car did not make finish line, no lap time!!')
        if isfield(daq,'lapInfo')
            daq = rmfield(daq,'lapInfo');
        end
    else    
        lapTime = diff(daq.rawData.time.meas([indStart indEnd]));
        daq.lapInfo.laps{1}.ind = [indStart indEnd];
        daq.lapInfo.laps{1}.lapTime = createDaqChannelData(lapTime,'s','Lap Time');    
    end
    
    
    %Let user know what's going on
    if ~suppressOutput
        disp(' ')
        fprintf('File: %s\n',displayDaqFiles(daq,'suppressOutput',true));
        fprintf('Starting distance = %5.3f at index %i\n',startDistance,indStart);
        fprintf('Finish distance = %5.3f at index %i\n',endDist,indEnd);
        fprintf('Calculated lap time = %5.3f [s]\n',lapTime);
        disp(' ')
    end
    
    %Update DAQ
    DAQ{iFile} = daq;
    
end

times = displayLapTimes(DAQ,'suppressOutput',suppressOutput);
