function daq = gpopsMPC(daq)
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
    masterDaq.status.currentDistance = masterDaq.header.s0;
    masterDaq.status.currentX0 = masterDaq.header.x0;
    masterDaq.status.currentSegment = 1;
    masterDaq.status.currentGuess = masterDaq.header.setup.guess;
  
else %Relaunch?
    error('need to code')
    masterDaq = daq;
    
    segDaq.header = masterDaq.header; 
end

%% Define useful variables - to clean up code
variableNames = masterDaq.header.variableNames;

%% Main Loop
checkeredFlag = false;
while ~checkeredFlag
    %% SETUP AND RUN
    
    %Setup the horizon
    segDaq = [];
    segDaq.header = masterDaq.header;   
    
    %Update the segDaq bounds for the ocp
    setup = segDaq.header.setup;
    setup.bounds.phase.initialtime.lower  = masterDaq.status.currentDistance; 
    setup.bounds.phase.initialtime.upper  = masterDaq.status.currentDistance;
    setup.bounds.phase.finaltime.lower    = masterDaq.status.currentDistance + masterDaq.header.horizon; 
    setup.bounds.phase.finaltime.upper    = masterDaq.status.currentDistance + masterDaq.header.horizon;
    setup.bounds.phase.initialstate.lower = masterDaq.status.currentX0 ;
    setup.bounds.phase.initialstate.upper = masterDaq.status.currentX0 ;
    setup.guess                           = masterDaq.status.currentGuess;
    segDaq.header.setup = setup;
    
    %Run the segment daq
    fprintf('HORIZON: %03i currently running....\n',masterDaq.status.currentSegment);
    [segDaq, convergence] = fourwheelMain(segDaq);
    
    %Save a snapshot of everything
    snapshotFilename = sprintf('%s_solutionSnapshot_Horizon%03i',datestr(now,'yyyy-mm-dd_HH_MM_SS'),masterDaq.status.currentSegment);
    save(snapshotFilename)
    
    %See if the problem converged
    if ~convergence
        warning('OCP failed to converge at horizon %03i',masterDaq.status.currentSegment);
        return
    end
    
    %% Update Master Solution
    %If it converged we got here. Next we need to update the master
    %solution. First grab the data just over the MPC update interval
    ind = findDaqIndForConditions(['distance>=' num2str(masterDaq.status.currentDistance) ...
                                   '& distance<=' num2str(masterDaq.status.currentDistance + masterDaq.header.controlHorizon)],...
                                   segDaq);
                               
    mpcIntervalDaq = assembleNewDaqAtIndicies(ind,segDaq,'normalizeIndepVarChannels',false);
    
    %Put it in the daq file
    if ~isfield(masterDaq,'rawData')                                       %If there isn't a rawData field, we need to start it
        masterDaq.rawData = mpcIntervalDaq.rawData;
    else
        channels = daqChannels(segDaq);
        for iCh = 1:size(channels,1) 
            ch = channels{iCh};
            masterDaq.rawData.(ch).meas = [masterDaq.rawData.(ch).meas; mpcIntervalDaq.rawData.(ch).meas];
        end
    end
    
    %Now grab the last index as the starting point of the next
    indNewX0 = length(masterDaq.rawData.(variableNames.indepVarName).meas);
    newX0daq = assembleNewDaqAtIndicies(indNewX0,masterDaq,'normalizeIndepVarChannels',false);
    newX0 = writeDaqChannelsToMatrix(newX0daq,'selectedChannels',variableNames.stateNames);
    
    %Now, remove the last point from the master daq as it'll be the first
    %point in the next solution
    masterDaq = removeDataFromDaqFileAtIndex(indNewX0,masterDaq);
    
    %Update the master solution
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
    
%       %Based off previous horizon:
%     nextGuess.phase.time    = segDaq.gpopsOutput.result.solution.phase.time;
%     nextGuess.phase.state   = segDaq.gpopsOutput.result.solution.phase.state;
%     nextGuess.phase.control = segDaq.gpopsOutput.result.solution.phase.control;
%     nextGuess.phase.integral= segDaq.gpopsOutput.result.solution.phase.integral;
    
        
      %Based off previous horizon:
    nextGuess.phase.time    = [segDaq.gpopsOutput.result.solution.phase.time; masterDaq.status.currentDistance + masterDaq.header.horizon];
    nextGuess.phase.state   = [segDaq.gpopsOutput.result.solution.phase.state; segDaq.gpopsOutput.result.solution.phase.state(end,:)];
    nextGuess.phase.control = [segDaq.gpopsOutput.result.solution.phase.control; 0];
    nextGuess.phase.integral= segDaq.gpopsOutput.result.solution.phase.integral;
    
    
%     %Arbitary initial guess:
%     nextGuess.phase.time    = [masterDaq.status.currentDistance; segDaq.rawData.distance.meas(end); masterDaq.status.currentDistance + masterDaq.header.horizon];
%     nextGuess.phase.state   = ones(3,1)*newX0;
%     nextGuess.phase.control = zeros(3,1);
%     nextGuess.phase.control(2) = 
%     nextGuess.phase.integral= 0;
    
    
    
    
    masterDaq.status.currentGuess = nextGuess;
   
    
    %See if we should end lap sim
    if masterDaq.status.currentDistance > masterDaq.header.finishDistance;
        checkeredFlag = true;
    end
end %While no checkered flag

% [times, masterDaq] = calculateManeuveringTime(masterDaq);
%Save the final solution
snapshotFilename = sprintf('%s_MasterDaqSolution_t=%05.3f',datestr(now,'yyyy-mm-dd_HH_MM_SS',times));
daq = masterDaq;
save(snapshotFilename,'daq');

    