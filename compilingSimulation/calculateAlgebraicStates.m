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

input.phase.time = output.result.(solutionField).phase.time;
input.phase.state = output.result.(solutionField).phase.state;
input.phase.control = output.result.(solutionField).phase.control;

input.auxdata = daq.header.setup.auxdata;

%Need to turn scaling off as the output is already unscaled
fields = fieldnames(input.auxdata.scaling);
for iField = 1:length(fields)
    field = fields{iField};
    input.auxdata.scaling.(field) = ones(size(input.auxdata.scaling.(field)));
end

phaseout = fourwheelContinuous(input);

algebraicStates = fieldnames(phaseout.algebraicStates);
for i = 1:length(algebraicStates)
    ch = algebraicStates{i};
    daq.rawData = addMathChannelsThatAreStandardChannels(daq.rawData,ch,phaseout.algebraicStates.(ch).meas,'','','suppressWarning',true);
end