function [horizonDaq, convergence] = fourwheelMain(horizonDaq)
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

setup = horizonDaq.header.setup;
setup.functions.continuous        = @fourwheelContinuous;
setup.functions.endpoint          = @fourwheelEndpoint;

adigatorfilenames = adigatorGenFiles4gpops2(setup);

setup.adigatorgrd.continuous    = @fourwheelContinuousADiGatorGrd;
setup.adigatorgrd.endpoint      = @fourwheelEndpointADiGatorGrd;
setup.adigatorhes.continuous    = @fourwheelContinuousADiGatorHes;
setup.adigatorhes.endpoint      = @fourwheelEndpointADiGatorHes;


%% SOVLE
tic
output   = gpops2(setup);
toc

saveSnapshotofShortSeg = 'snapshot';
%%Deal w/ saving file
if ~isempty(saveSnapshotofShortSeg)
    save(saveSnapshotofShortSeg)
    ipoptfn = sprintf('%s_IPOPT.txt',strrep(saveSnapshotofShortSeg,'.mat',''));
    sysCommand = sprintf('mv quadCarIPOPTinfo.txt %s',ipoptfn);
    system(sysCommand);
end

%Post process
if ismember(output.result.nlpinfo,horizonDaq.header.acceptableNlpOutpus)
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
horizonDaq = calculateAlgebraicStates(horizonDaq);


%%Deal w/ saving file - need to resave after editing
if ~isempty(saveSnapshotofShortSeg)
    save(saveSnapshotofShortSeg)
end

        

