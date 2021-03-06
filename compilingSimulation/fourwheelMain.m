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
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

if ~exist('horizonDaq','var')
    horizonDaq = generateInitialDaq();
end


defaults = {'saveSnapshotofShortSeg','snapshot'
            'calcAlgebraicStates',true};
setDefaultsForVarargin(defaults,varargin);

setup = horizonDaq.header.setup;

setup.functions.continuous        = @fourwheelContinuous;                  %Unfortunately, I don't know a better way than a static link to the continuous and endpoitn functions.
setup.functions.endpoint          = @fourwheelEndpoint;
% adigatorfilenames = adigatorGenFiles4gpops2(setup);
setup.adigatorgrd.continuous    = @fourwheelContinuousADiGatorGrd;         %Unfortunately, I don't know a better way than a static link to the adigator files.
setup.adigatorgrd.endpoint      = @fourwheelEndpointADiGatorGrd;
setup.adigatorhes.continuous    = @fourwheelContinuousADiGatorHes;
setup.adigatorhes.endpoint      = @fourwheelEndpointADiGatorHes;


%% SOVLE
tic
output   = gpops2(setup);
toc


%%Deal w/ saving file
if ~isempty(saveSnapshotofShortSeg)
    save(saveSnapshotofShortSeg)
    ipoptfn = sprintf('%s_IPOPT.txt',strrep(saveSnapshotofShortSeg,'.mat',''));
    sysCommand = sprintf('mv quadCarIPOPTinfo.txt %s',ipoptfn);
    system(sysCommand);
end

%Post process
if ismember(output.result.nlpinfo,horizonDaq.header.acceptableNlpOutputs)
    convergence = true;
%     if ~isa(gpopsOptions.specifyMeshIterationSolution,'char')
%         output.result = output.meshhistory(gpopsOptions.specifyMeshIterationSolution).result;
%     end
else
    horizonDaq = [];
    convergence = false; 
    return
end

%horizonDaq the solution
                            %    vx   vy   r  t
daqGpops = parseGpops2toDaq(output,...
    horizonDaq.header.variableNames.stateNames,...
    horizonDaq.header.variableNames.controlNames,...
    horizonDaq.header.variableNames.indepVarName,...
    horizonDaq.header.variableNames.units,...
    horizonDaq.header.variableNames.names);
tempHeader = horizonDaq.header;                          
horizonDaq = catstruct(daqGpops,horizonDaq);
horizonDaq.gpopsSetup = setup;
horizonDaq.header = tempHeader;

%Calc algebraic states
if calcAlgebraicStates
    horizonDaq = calculateAlgebraicStates(horizonDaq);
end
horizonDaq.rawData.(horizonDaq.header.variableNames.indepVarName).meas = ...
    horizonDaq.rawData.(horizonDaq.header.variableNames.indepVarName).meas + setup.auxdata.currentDistance;


%%Deal w/ saving file - need to resave after editing
if ~isempty(saveSnapshotofShortSeg)
    save(saveSnapshotofShortSeg)
end

        

