function daq = gpopsMPC(daq)
%This funciton performs the MPC lap simulation with each horizon using
%GPOPS-II to solve the OCP over the horizon.
%INTPUS: 
%    daq - see generateInitialDaq();
%OUTPUS:
%    daq - updated daq file
%Creation: 21 Dec 2017 - Jeff Anderson 
%Updated:  02 May 2018 - Jeff Anderson - Added a return if sim is already 
%    done to help out batch running.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%Fresh sim?
if nargin == 0 
    masterDaq = generateInitialDaq();    
else
    masterDaq = daq;
end

%Make sure the sim hasn't finished already
if masterDaq.header.simFinished
    fprintf('Simulation already complete\n')
    return
end

if ~isfield(masterDaq,'status')
    masterDaq.status.currentDistance = masterDaq.header.initialDistance;
    masterDaq.status.currentX0 = masterDaq.header.x0;
    masterDaq.status.currentSegment = 1;
    masterDaq.status.currentGuess = masterDaq.header.setup.guess;
end

%% Define useful variables - to clean up code
variableNames = masterDaq.header.variableNames;

%% Main Loop
checkeredFlag = false;
while ~checkeredFlag

    genStats(masterDaq);
    %% Horizon refinement
    segDaq = [];

    convergence = false;
    iHorizonRefinement = 1;
    segDaq.header.iHorizonRefinement = iHorizonRefinement;
    horizon = masterDaq.header.horizon;
    x0 = masterDaq.status.currentX0;
    
    %Setup diary
    if masterDaq.header.saveDiaryFiles
        diaryFilename = sprintf('Horizon%03i-Diary',masterDaq.status.currentSegment);
        diary([diaryFilename '_diary']);
    end
    
    %Run the sim until converged
    while ~convergence
        %Update the segDaq bounds for the ocp
        segDaq.header = masterDaq.header; 
        setup = segDaq.header.setup;
        setup.bounds.phase.initialtime.lower  = 0*setup.auxdata.scaling.length;
        setup.bounds.phase.initialtime.upper  = 0*setup.auxdata.scaling.length;
        setup.bounds.phase.finaltime.lower    = horizon*setup.auxdata.scaling.length;
        setup.bounds.phase.finaltime.upper    = horizon*setup.auxdata.scaling.length;
        setup.bounds.phase.initialstate.lower = x0.*setup.auxdata.scaling.state ;
        setup.bounds.phase.initialstate.upper = x0.*setup.auxdata.scaling.state ;
        setup.guess                           = masterDaq.status.currentGuess;

        %Auxdata to get current distance
        setup.auxdata.currentDistance         = masterDaq.status.currentDistance;
        setup.auxdata.costTime                = masterDaq.header.switchingDaq.rawData.switching.meas(masterDaq.status.currentSegment,1);
        setup.auxdata.costVx                  = masterDaq.header.switchingDaq.rawData.switching.meas(masterDaq.status.currentSegment,2);
        setup.auxdata.vehicle                 = masterDaq.vehicle;
        segDaq.header.setup = setup;        
        
        fprintf('HORIZON: %03i currently running....\n',masterDaq.status.currentSegment);
      
        %Run the segment daq
        
        [segDaq, convergence] = fourwheelMain(segDaq,'calcAlgebraicStates',false,'saveSnapshotofShortSeg',[]);
        
        %If horizon refinement is on, decrease horizon here
        if ~convergence && ...
           masterDaq.header.horizonRefinement && ...
           horizon-masterDaq.header.horizonDecrement >= masterDaq.header.minimumHorizon
       
            %Increment counters
            iHorizonRefinement = iHorizonRefinement+1;
            segDaq.header.iHorizonRefinement = iHorizonRefinement;
            
            %Shorten the horizon
            horizon = horizon - masterDaq.header.horizonDecrement;
            segDaq.header.horizon = horizon;
            
        %If the simulation failed and we've already refined as much as
        %possible, kick out
        elseif ~convergence 
            masterDaq.header.simFinished = true;
            masterDaq.header.conv = false;
            genStats(masterDaq);
            fprintf('Simulation failed on Horizon %03i',masterDaq.status.currentSegment)
            return
        end
    end%while ~converge
    
    %Save a snapshot of everything
    fprintf('HORIZON: %03i EXIT, convergence = %d.\n',masterDaq.status.currentSegment,convergence);
    
    %Fix up diary file
    if masterDaq.header.saveDiaryFiles
        diary off
    end
    
    %Save the OCP over the horizon, first clean out the header (files
    %getting big)
    justHorizonFilename = sprintf('Horizon%03i-OCP',masterDaq.status.currentSegment);
    tempDaq = segDaq;
    daq = tempDaq;
    daq = rmfield(daq,'header');
    daq.header.filename = tempDaq.header.filename;
    daq.header.path = tempDaq.header.path;
    save(justHorizonFilename,'daq');
    
   
    %% Update Master Solution
    %If it converged we got here. Next we need to update the master
    %solution. First grab the data just over the MPC update interval
    ind = [1:findNearestPoint(masterDaq.status.currentDistance + masterDaq.header.controlHorizon,segDaq.rawData.distance.meas)];
    mpcIntervalDaq = assembleNewDaqAtIndicies(ind,segDaq,'normalizeIndepVarChannels',false);

    %Put it in the daq file
    if ~isfield(masterDaq,'rawData')                                       %If there isn't a rawData field, we need to start it
        masterDaq.rawData = mpcIntervalDaq.rawData;
        masterDaq.rawData = addMathChannelsThatAreStandardChannels(masterDaq.rawData,'mpcHorizon',ones(size(mpcIntervalDaq.rawData.(variableNames.indepVarName).meas)));
    else
        channels = daqChannels(segDaq);
        for iCh = 1:size(channels,1) 
            ch = channels{iCh};
            masterDaq.rawData.(ch).meas = [masterDaq.rawData.(ch).meas; mpcIntervalDaq.rawData.(ch).meas];
        end
        masterDaq.rawData.mpcHorizon.meas =  [masterDaq.rawData.mpcHorizon.meas; masterDaq.status.currentSegment*ones(size(mpcIntervalDaq.rawData.(variableNames.indepVarName).meas))];
    end
    
    
    %Now grab the last index as the starting point of the next
    indNewX0 = length(masterDaq.rawData.(variableNames.indepVarName).meas);
    newX0daq = assembleNewDaqAtIndicies(indNewX0,masterDaq,'normalizeIndepVarChannels',false);
    newX0 = writeDaqChannelsToMatrix(newX0daq,'selectedChannels',variableNames.stateNames);
    
    %Now, remove the last point from the master daq as it'll be the first
    %point in the next solution
    masterDaq = removeDataFromDaqFileAtIndex(indNewX0,masterDaq);
    
    %Save Update the master solution 
    masterDaq.status.currentDistance = masterDaq.status.currentDistance + masterDaq.header.controlHorizon;
    masterDaq.status.currentX0 = newX0;
    masterDaq.status.currentSegment = masterDaq.status.currentSegment + 1;
    
    %% Next GUESS
    %Need to establish the next horizon's initial guess based on the last
    %horizon. Start by zeros and then import the last horizon over it. So
    %we have the correct number of trailing zeros
%     sGuess = (masterDaq.status.currentDistance:...
%               masterDaq.header.interpolationAccuracy:...
%               masterDaq.status.currentDistance + masterDaq.header.horizon)';
%     uGuess = zeros(length(sGuess),length(variableNames.controlNames)); %pre allocate
    
    %Grab last seg controls and interp to the guess
%     uLastSeg = writeDaqChannelsToMatrix(segDaq,'selectedChannels',variableNames.controlNames);
%     distLastSeg = writeDaqChannelsToMatrix(segDaq,'selectedChannels',variableNames.indepVarName);
%     for iCh = 1:length(masterDaq.header.variableNames.controlNames)
%         uGuess(:,iCh) = interp1(distLastSeg,uLastSeg,sGuess,'spline',0);
%     end


    
    %Propogate the dynamics and establish the next guess
%     guessDaq.header.setup = masterDaq.header.setup; %Mimic how we'll have the daq in the final form for OCP
%     guessDaq = generateGuessDaq(sGuess,newX0,uGuess,guessDaq);
%     nextGuess.phase.time    = writeDaqChannelsToMatrix(guessDaq,'selectedChannels',variableNames.indepVarName);
%     nextGuess.phase.state   = writeDaqChannelsToMatrix(guessDaq,'selectedChannels',variableNames.stateNames);
%     nextGuess.phase.control = writeDaqChannelsToMatrix(guessDaq,'selectedChannels',variableNames.controlNames);
%     nextGuess.phase.integral     = 0;    
    
        
%       %Based off previous horizon: grab the last horizon, shift for the mpc
%       %update interval and then tack on some zeros
%       indOfLastHorizonForGuess = findNearestPoint(segDaq.gpopsOutput.result.solution.phase.time,masterDaq.header.controlHorizon);
%       nextGuess.phase.time    = [segDaq.gpopsOutput.result.solution.phase.time(indOfLastHorizonForGuess:end)-segDaq.gpopsOutput.result.solution.phase.time(indOfLastHorizonForGuess);  masterDaq.header.horizon];
%       nextGuess.phase.state   = [segDaq.gpopsOutput.result.solution.phase.state(indOfLastHorizonForGuess:end,:); segDaq.gpopsOutput.result.solution.phase.state(end,:)];
%       nextGuess.phase.control = [segDaq.gpopsOutput.result.solution.phase.control(indOfLastHorizonForGuess:end,:); zeros(1,length(masterDaq.header.variableNames.controlNames))];
%       nextGuess.phase.integral= segDaq.gpopsOutput.result.solution.phase.integral;
      
      %Make it a sparse guess
%       nextGuess.phase.time    = [segDaq.gpopsOutput.result.solution.phase.time(indOfLastHorizonForGuess)-segDaq.gpopsOutput.result.solution.phase.time(indOfLastHorizonForGuess); segDaq.gpopsOutput.result.solution.phase.time(end)-segDaq.gpopsOutput.result.solution.phase.time(indOfLastHorizonForGuess);  masterDaq.header.horizon];
%       nextGuess.phase.state   = [segDaq.gpopsOutput.result.solution.phase.state(indOfLastHorizonForGuess,:); segDaq.gpopsOutput.result.solution.phase.state(end,:); segDaq.gpopsOutput.result.solution.phase.state(end,:)];
%       nextGuess.phase.control = [segDaq.gpopsOutput.result.solution.phase.control(indOfLastHorizonForGuess,:); segDaq.gpopsOutput.result.solution.phase.control(end,:); zeros(1,length(masterDaq.header.variableNames.controlNames))];
      

%     %Arbitary initial guess:
%     nextGuess.phase.time    = [masterDaq.status.currentDistance; segDaq.rawData.distance.meas(end); masterDaq.status.currentDistance + masterDaq.header.horizon];
%     nextGuess.phase.state   = ones(3,1)*newX0;
%     nextGuess.phase.control = zeros(3,1);
%     nextGuess.phase.control(2) = 
%     nextGuess.phase.integral= 0;
    
    
    
    
%     masterDaq.status.currentGuess = nextGuess;
   
    
    %See if we should end lap sim
    if masterDaq.status.currentDistance > masterDaq.header.finishDistance;
        checkeredFlag = true;
    end
    
    %Save the master daq
    daq = masterDaq;
    save(masterDaq.header.filename,'daq');
    
    
   
end %While no checkered flag

masterDaq.header.simFinished = true;
masterDaq.header.conv = true;
[~, masterDaq] = calculateManeuveringTime(masterDaq);

genStats(masterDaq);

%Save the final solution
daq = masterDaq;
save(masterDaq.header.filename,'daq');


    