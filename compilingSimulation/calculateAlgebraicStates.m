function daq = calculateAlgebraicStates(daq)
%This funciton will throw the inputs and states back through the ODE and
%calculate the algebraic states in the fourwheelContinuous model from
%Gpops2.
%INPTUS:
%    daq (optional) - local instance of the daq file
%OUTPUTS:
%    daq (optional) - updated local instnace of the daq file
%
%Creation: 12 Nov 2017 - Jeff Anderson
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


output = daq.gpopsOutput;

solutionField = 'interpsolution'; %field in the gpops results output.results.(solutionField) originally just solution but looks like the interpsol is more what I need.

input.phase.time = output.result.(solutionField).phase.time - output.result.(solutionField).phase.time(1);
input.phase.state = output.result.(solutionField).phase.state;
input.phase.control = output.result.(solutionField).phase.control;

input.auxdata = daq.header.setup.auxdata;

%Need to turn scaling off as the output is already unscaled
fields = fieldnames(input.auxdata.scaling);
for iField = 1:length(fields)
    field = fields{iField};
    input.auxdata.scaling.(field) = ones(size(input.auxdata.scaling.(field)));
end

%% CAll CONTINOUS FUNCTION
%combine the continuous function and algebraic state command files into one
%mcode
if ~isunix; error('Following code uses bash, need to recode for windows'); end
!rm temp.m
!echo 'function phaseout = temp(input)' > temp.m
!tail -n +2 fourwheelContinuous.m >> temp.m
!cat algebraicStatesCommandsForContinuousFunction.m >> temp.m

%Sometimes on palmetto, writting files takes time...so need try more 1x
count = 1; 
retryFlag = true;
err = [];
while count<100 && retryFlag
    try
        phaseout = temp(input);
    catch err
    end
    if isempty(err)
        retryFlag = false;
    else
        count = count+1;
        err = [];
    end
end
!rm temp.m
algebraicStates = fieldnames(phaseout.algebraicStates);
S = daq.rawData.distance.meas-daq.rawData.distance.meas(1);
for i = 1:length(algebraicStates)
    ch = algebraicStates{i};
    
    %Need to reinterpolate them to the output
    data = interp1(input.phase.time,phaseout.algebraicStates.(ch).meas,S);
    daq.rawData = addMathChannelsThatAreStandardChannels(daq.rawData,ch,data,'','','suppressWarning',true);
end