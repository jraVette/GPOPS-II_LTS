function [horizonDaq, convergence] = fourwheelMain(horizonDaq,varargin)
%This funciton will run a short segment OCP
%
%INPTUS:
%   daq - structure of the horizon with all necessary GPOPS fields.
%   convergence - true/false whether the segment converged accourding to
%       NLP output and acceptable solutions defined in 
%       daq.header.acceptableNlpOutpus.
%
%OUTPUTS:
%   daq - updated strucutre with the solution
%
%Creation: 21 Dec 2017 - Jeff Anderson
%Updated:  25 Feb 2018 - Jeff Anderson - have it returing the failed output
%   vs an empty horizonDaq. Killed the snapshot code as the diary is
%   handled in gpopsMPC.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
defaults = {'retryIfOneMeshSucessful',true
            'nRetries', 2};
setDefaultsForVarargin(defaults,varargin);

setup = horizonDaq.header.setup;

setup.functions.continuous        = @fourwheelContinuous;                  %Unfortunately, I don't know a better way than a static link to the continuous and endpoitn functions.
setup.functions.endpoint          = @fourwheelEndpoint;

setup.adigatorgrd.continuous    = @fourwheelContinuousADiGatorGrd;         %Unfortunately, I don't know a better way than a static link to the adigator files.
setup.adigatorgrd.endpoint      = @fourwheelEndpointADiGatorGrd;
setup.adigatorhes.continuous    = @fourwheelContinuousADiGatorHes;
setup.adigatorhes.endpoint      = @fourwheelEndpointADiGatorHes;


%% SOVLE
output   = gpops2(setup);

nlpInfoOverMeshHistories = zeros(output.meshcounts,1);
for iMesh = 1:output.meshcounts
    nlpInfoOverMeshHistories(iMesh) = output.meshhistory(iMesh).result.nlpinfo;
end

%If it didn't converge, see if any mesh histories did, and if so try them
%as a guess.
retryCount =1 ;
while ~ismember(output.result.nlpinfo,horizonDaq.header.acceptableNlpOutpus) && ...
      ~isempty(find(nlpInfoOverMeshHistories == 0, 1)) && ...
      retryIfOneMeshSucessful && ...
      retryCount <= nRetries 
    %Didn't converge, but we had one mesh history that did, use that as a
    %guess
    iMeshForGuess = find(nlpInfoOverMeshHistories == 0,1,'last');
    setup.guess.phase.time = output.meshhistory(iMeshForGuess).result.solution.phase.time;
    setup.guess.phase.state = output.meshhistory(iMeshForGuess).result.solution.phase.state;
    setup.guess.phase.control = output.meshhistory(iMeshForGuess).result.solution.phase.control;
    setup.guess.phase.integral = output.meshhistory(iMeshForGuess).result.solution.phase.integral;
    output   = gpops2(setup);
    retryCount = retryCount+1;
end


%Post process
if ismember(output.result.nlpinfo,horizonDaq.header.acceptableNlpOutpus)
    convergence = true;
    %horizonDaq the solution
                                %    vx   vy   r  t
    daqGpops = parseGpops2toDaq(output,...
        horizonDaq.header.variableNames.stateNames,...
        horizonDaq.header.variableNames.controlNames,...
        horizonDaq.header.variableNames.indepVarName,...
        horizonDaq.header.variableNames.units,...
        horizonDaq.header.variableNames.names,...
        'dIndependantVariable',horizonDaq.header.interpolationAccuracy);
    tempHeader = horizonDaq.header;                          
    horizonDaq = catstruct(daqGpops,horizonDaq);
    horizonDaq.gpopsSetup = setup;
    horizonDaq.header = tempHeader;

    %Calc algebraic states
    horizonDaq = calculateAlgebraicStates(horizonDaq);
    
    %Update the distance to the actual distance
    horizonDaq.rawData.distance.meas = horizonDaq.rawData.distance.meas + setup.auxdata.currentDistance;    
else
    horizonDaq.gpopsOutput = output;
    horizonDaq.gpopsSetup = setup;
    convergence = false; 
end


        

