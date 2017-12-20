setup = daq.header.setup;
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
if ismember(output.result.nlpinfo,daq.header.acceptableNlpOutpus)
    convergence = true;
%     if ~isa(gpopsOptions.specifyMeshIterationSolution,'char')
%         output.result = output.meshhistory(gpopsOptions.specifyMeshIterationSolution).result;
%     end
else
    daq = [];
    convergence = false; 
    return
end

%Daq the solution
                            %    vx   vy   r  t
daqGpops = parseGpops2toDaq(output,daq.header.variableNames.stateNames,...
                              daq.header.variableNames.controlNames,...
                              daq.header.variableNames.indepVarName,...
                              daq.header.variableNames.units,...
                              daq.header.variableNames.names);
tempHeader = daq.header;                          
daq = catstruct(daqGpops,daq);
daq.gpopsSetup = setup;
daq.header = tempHeader;

%Calc algebraic states
daq = calculateAlgebraicStates(daq);


%%Deal w/ saving file - need to resave after editing
if ~isempty(saveSnapshotofShortSeg)
    save(saveSnapshotofShortSeg)
end

        

