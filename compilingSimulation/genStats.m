function genStats(daq)
%This will generate stats on the current simulations
%Updated: 23 Jan 2017 - Jeff Anderson - I added sim finished to the daq
%file in hopes of catching non convergent runs.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


%% Re load from existing file or start a new one
if exist('stat.mat','file')
    load('stat.mat')
else
    stat.simStarted         = now;
    stat.simStartedString   = datestr(stat.simStarted);
    stat.lastCheckin        = now;
    stat.lastCheckinString  = datestr(stat.lastCheckin);
    stat.elapsedTime        = now-now;
    stat.elapsedTimeString  = datestr(stat.elapsedTime,'HH:MM:SS');
    stat.currentHorizon     = 1;
    stat.simFinished        = 0;
    stat.lapTime            = nan;
    stat.currentDistance    = 0;
    stat.conv               = nan;
    stat.iterNumb           = nan;
    stat.filename           = displayDaqFiles(daq,'suppressOutput',true,'lookForShortFilename',false);
    stat.directory          = daq.header.path;
end


%% Deal with the daq file

if exist('daq','var') 
    %Current horizon
    stat.currentHorizon = daq.status.currentSegment;
    
    %Update the distance
    distance = [];
    getChannelDataFromDaqFile(daq,'distance');
    if ~isempty(distance)
        stat.currentDistance = distance(end);
    end
    
    %See if simulation is finished
    if isfield(daq.header,'simFinished')
        stat.simFinished = daq.header.simFinished;
    end
    
    %Update the convergence
    if isfield(daq.header,'conv')
        stat.conv = daq.header.conv;
    end
    
    %Update laptime
    times = displayLapTimes(daq);
    if ~isempty(times)
        stat.lapTime = times;
    end
end

%% Deal with time since start:
stat.elapsedTime       = now-stat.lastCheckin+stat.elapsedTime;
stat.lastCheckin       = now;
stat.lastCheckinString = datestr(stat.lastCheckin);
stat.elapsedTimeString = datestr(stat.elapsedTime,'HH:MM:SS');
save('stat.mat','stat')
fprintf('Elapsed time: %s\n',stat.elapsedTimeString);




%% Setup a text file
statFile = 'stat.txt';
if exist(statFile,'file')
    system(sprintf('rm %s',statFile));
end

fileID = fopen(statFile,'w');
fprintf(fileID,'Iterate             = %i\n',   stat.iterNumb);
fprintf(fileID,'Sim Started         = %s\n',   stat.simStartedString);
fprintf(fileID,'Last Check In       = %s\n',   stat.lastCheckinString);
fprintf(fileID,'Elapsed Time        = %s\n',   stat.elapsedTimeString);
fprintf(fileID,'Current Distance    = %5.1f\n',stat.currentDistance);
fprintf(fileID,'Simulation Finished = %1i\n',  stat.simFinished);
fprintf(fileID,'Lap Time            = %6.4f\n',  stat.lapTime);
fprintf(fileID,'Convergence         = %f\n',     stat.conv);
fclose(fileID);
    
    
