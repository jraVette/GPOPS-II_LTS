%-------------------------------------------%
% BEGIN: function particleMotionEndpoint.m %
%-------------------------------------------%
function output = fourwheelEndpoint(input)
% cost = input.auxdata.cost;
% %Objective
integral = input.phase(1).integral;
% terminalCosts = sum(((input.auxdata.terminalCost.targetState - input.phase.finalstate)./input.auxdata.terminalCost.scaling).^2)/length(input.phase.finalstate) ;
% output.objective = integral + terminalCosts*input.auxdata.terminalCost.weight;
output.objective = integral;

end
%-------------------------------------------%
% END: function brachistochroneEndpoint.m   %
%-------------------------------------------%
