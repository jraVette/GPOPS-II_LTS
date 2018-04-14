function daq = gpopsMPC_horizonRefinement3(daq)
%This funciton performs the MPC lap simulation with each horizon using
%GPOPS-II to solve the OCP over the horizon.
%INTPUS: 
%    daq - see generateInitialDaq();
%OUTPUS:
%    daq - updated daq file
%Creation: 21 Dec 2017 - Jeff Anderson
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%Fresh sim?
if nargin == 0 
    masterDaq = generateInitialDaq();    
else
    masterDaq = daq;
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
    diaryFilename = sprintf('Horizon%03i-Diary',masterDaq.status.currentSegment);
    diary([diaryFilename '_diary']);
    
    while ~convergence
        %Update the segDaq bounds for the ocp
        segDaq.header = masterDaq.header; 
        setup = segDaq.header.setup;
        setup.bounds.phase.initialtime.lower  = 0*setup.auxdata.scaling.length;1
        setup.bounds.phase.initialtime.upper  = 0*setup.auxdata.scaling.length;
        setup.bounds.phase.finaltime.lower    = horizon*setup.auxdata.scaling.length;
        setup.bounds.phase.finaltime.upper    = horizon*setup.auxdata.scaling.length;
        setup.bounds.phase.initialstate.lower = x0.*setup.auxdata.scaling.state ;
        setup.bounds.phase.initialstate.upper = x0.*setup.auxdata.scaling.state ;
        setup.guess                           = masterDaq.status.currentGuess;

        %Auxdata to get current distance
        setup.auxdata.currentDistance         = masterDaq.status.currentDistance;
        segDaq.header.setup = setup;        
        
        fprintf('HORIZON: %03i currently running....\n',masterDaq.status.currentSegment);
      
        %Run the segment daq
        
        [segDaq, convergence] = fourwheelMain(segDaq,'calcAlgebraicStates',false,'saveSnapshotofShortSeg',[]);
        
        if ~convergence && iHorizonRefinement <= 3 && masterDaq.status.currentSegment > 1
            iHorizonRefinement = iHorizonRefinement+1;
            segDaq.header.iHorizonRefinement = iHorizonRefinement;
            
            %Back the solution up and lengthen the horizon
            masterDaq.status.currentDistance = masterDaq.status.currentDistance - masterDaq.header.controlHorizon/4;
%             horizon = horizon + masterDaq.header.controlHorizon;
            
            %Find the new x0
            indCurrentDistance = find(masterDaq.status.currentDistance == masterDaq.rawData.distance.meas);
            newX0daq = assembleNewDaqAtIndicies(indCurrentDistance,masterDaq,'normalizeIndepVarChannels',false);
            x0 = writeDaqChannelsToMatrix(newX0daq,'selectedChannels',variableNames.stateNames);
            masterDaq.status.currentX0 = x0;
            
            %Rewind the master daq solution
            masterDaq = removeDataFromDaqFileAtIndex(indCurrentDistance:length(masterDaq.rawData.distance.meas),masterDaq);
            
            %Update the guess
            previousHorizonDaq = load(sprintf('Horizon%03i-OCP.mat',masterDaq.status.currentSegment-1));
            previousHorizonDaq = previousHorizonDaq.daq;
            lastHorizonDistances = previousHorizonDaq.gpopsOutput.result.solution.phase.time+previousHorizonDaq.gpopsSetup.auxdata.currentDistance;
            indOfLastHorizonForGuess = findNearestPoint(lastHorizonDistances,masterDaq.status.currentDistance);
            nextGuess.phase.time    = [previousHorizonDaq.gpopsOutput.result.solution.phase.time(indOfLastHorizonForGuess:end)-previousHorizonDaq.gpopsOutput.result.solution.phase.time(indOfLastHorizonForGuess);  masterDaq.header.horizon];
            nextGuess.phase.state   = [previousHorizonDaq.gpopsOutput.result.solution.phase.state(indOfLastHorizonForGuess:end,:); previousHorizonDaq.gpopsOutput.result.solution.phase.state(end,:)];
            nextGuess.phase.control = [previousHorizonDaq.gpopsOutput.result.solution.phase.control(indOfLastHorizonForGuess:end,:); zeros(1,length(masterDaq.header.variableNames.controlNames))];
            
            nextGuess.phase.integral= previousHorizonDaq.gpopsOutput.result.solution.phase.integral;            
            masterDaq.status.currentGuess = nextGuess;
        elseif ~convergence && horizon >= masterDaq.header.horizonRefinementMinHorizon;
            iHorizonRefinement = iHorizonRefinement+1;
            segDaq.header.iHorizonRefinement = iHorizonRefinement;
            
            %Back the solution up and lengthen the horizon
            horizon = horizon - masterDaq.header.horizonRefinementDecrement;            
        elseif ~convergence
            warning('HORIZON REFINEMENT FAILED - EXITING')
            break
        end
        
        
    end
    
    %Save a snapshot of everything
    fprintf('HORIZON: %03i EXIT, convergence = %d.\n',masterDaq.status.currentSegment,convergence);
    diary off
%     snapshotFilename = sprintf('Horizon%03i-Snapshot',masterDaq.status.currentSegment);
%     save(snapshotFilename)
    justHorizonFilename = sprintf('Horizon%03i-OCP',masterDaq.status.currentSegment);
    daq = segDaq;
    save(justHorizonFilename,'daq');
    
    
    %See if the problem converged
    if ~convergence
        warning('OCP failed to converge at horizon %03i',masterDaq.status.currentSegment);
        return
    end
    
    %% Update Master Solution
    %If it converged we got here. Next we need to update the master
    %solution. First grab the data just over the MPC update interval
%     ind = findDaqIndForConditions(['distance>=' sprintf('%f',masterDaq.status.currentDistance) ...
%                                    '& distance<=' sprintf('%f',masterDaq.status.currentDistance + masterDaq.header.controlHorizon)],...
%                                    segDaq);
    ind = [1:findNearestPoint(masterDaq.status.currentDistance + masterDaq.header.controlHorizon,segDaq.rawData.distance.meas)];
    mpcIntervalDaq = assembleNewDaqAtIndicies(ind,segDaq,'normalizeIndepVarChannels',false);
%     mpcIntervalDaq.rawData.distance.meas = mpcIntervalDaq.rawData.distance.meas + masterDaq.status.currentDistance;
    
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
    
    %See when we skip a beat
    if length(unique(diff(masterDaq.rawData.distance.meas))) >1
        save('snap.mat')
        error('got some issues')
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
    masterDaqFilename = sprintf('Horizon%03i-MasterDaq',masterDaq.status.currentSegment-1);
    save(masterDaqFilename,'daq');
    
   
end %While no checkered flag


genStats(masterDaq);
% [times, masterDaq] = calculateManeuveringTime(masterDaq);
%Save the final solution
snapshotFilename = sprintf('%s_MasterDaqSolution_t',datestr(now,'yyyy-mm-dd_HH_MM_SS'));
daq = masterDaq;
save(snapshotFilename,'daq');

    