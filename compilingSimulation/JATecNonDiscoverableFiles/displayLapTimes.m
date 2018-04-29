function [allTimes, allLapTimesPerFile] = displayLapTimes(daq,varargin)
%This function will display the laptimes of the global DAQ structures.
%Made for a quick display
%INPTUS: 
%    daq - (optional) local instance of a daq file use [] for global files
%    varargin - used to set defaults varaible below
%OUTPUS:
%    allTimes - double array of all lap times recorded
%    allLapTimesPerFile - cell array of lap times per file
%
%Creation: 2015 Jan 19 - Jeff Anderson
%Updated: 14 Dec 2015 - Jeff Anderson - updated fcn w/ new jatec utilites
%    and made ouptut.
%Updated: 28 Jun 2017 - Jeff Anderson - updated laps to print with date
%    strings.
%Updated: 07 Aug 2017 - Jeff Anderson - fixed that errored if no laps were
%    found.
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

if ~exist('daq','var'); daq = []; end
DAQ = useLocalOrGlobalDaqFiles(daq);

defaults = {'suppressOutput',false};
setDefaultsForVarargin(defaults,varargin);

if ~suppressOutput
    disp('Lap times per file:')
    disp('___________________________________________________________________')
    disp(' ')
end

count = 1;
allTimes = []; %Incase there are none, need variable defined
allLapTimesPerFile = {}; %Incase there are none, need variable defined

for iFile = 1:length(DAQ)
    daq = DAQ{iFile};
    
    if isfield(daq,'lapInfo')                                              %Is there any lap info
        if ~suppressOutput
            fprintf('Laps for file: %s\n',daq.header.filename)                 %If so, tell the user what file
        end
        
        for iLap = 1:length(daq.lapInfo.laps);                            %Loop over laps
            if isfield(daq.lapInfo.laps{iLap},'lapTime')                   %Make sure there is a lap time field
                lapTime = daq.lapInfo.laps{iLap}.lapTime;                  %Grab daqMeasure object
                
                if ~suppressOutput
%                     fprintf('Lap %i : %6.3f [%s]\n',iLap,lapTime.meas,lapTime.units)
                    dateNumberOfLapTime = datenum(0,0,0,0,0,lapTime.meas);
                    fprintf('Lap %i : %s\n',iLap,datestr(dateNumberOfLapTime,'MM:SS.FFF'))
                end
                
                if isempty(lapTime.meas)
                    lapTime.meas = nan;
                end
                allLapTimesPerFile{iFile}(iLap) = lapTime.meas;
                
                
                %Record all the times and what lap and what file
                allTimes(count) = lapTime.meas;                            %#ok dyn growth, don't known apriori
                allTimeLapCount(count) = iLap;                             %#ok dyn growth, don't known apriori
                allTimesFile{count} = daq.header.filename;                 %#ok dyn growth, don't known apriori
                allTimesPath{count} = daq.header.path;
                allTimesDaqInd(count) = iFile;
                count = count+1;
            else
                if ~suppressOutput
                    fprintf('No lap time metric found for lap %i\n',iLap)
                end
                allLapTimesPerFile{iFile} = [];
            end
        end
        
    else
        allLapTimesPerFile{iFile} = [];
        if ~suppressOutput
            fprintf('No lap info found on file: %s\n',daq.header.filename)
        end
    end
    
    if ~suppressOutput
        disp(' ')
    end
end

%%Display min lap time
if ~isempty(allTimes)
    minTime = min(allTimes);
    ind = find(minTime == allTimes);
    
    if ~suppressOutput
        disp('Minimum laptime:')
        for i = 1:length(ind)
            disp('___________________________________________________________________')
%             fprintf('%6.3f [s], Lap %2i, File: %s\n',allTimes(ind(i)),allTimeLapCount(ind(i)),allTimesFile{ind(i)});
            dateNumberOfLapTime = datenum(0,0,0,0,0,allTimes(ind(i)));
            fprintf('*%9s* Lap %2i, File: %s\n',datestr(dateNumberOfLapTime,'MM:SS.FFF'),allTimeLapCount(ind(i)),allTimesFile{ind(i)});
            fprintf('            DAQ index: %i\n',allTimesDaqInd(ind(i)));
            fprintf('            Path: %s\n',allTimesPath{ind(i)});
            disp(' ')
            disp('Time [s]:')
        end
    end
else
    if isempty(DAQ)
        disp('No files loaded')
    else
        disp('No laps found')
    end
end 


