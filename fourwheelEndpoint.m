%-------------------------------------------%
% BEGIN: function particleMotionEndpoint.m %
%-------------------------------------------%
function output = fourwheelEndpoint(input)
% cost = input.auxdata.cost;
% %Objective

output.objective = -input.phase(1).finalstate(1) + input.auxdata.controlWeight*input.phase(1).integral;

end
%-------------------------------------------%
% END: function brachistochroneEndpoint.m   %
%-------------------------------------------%
