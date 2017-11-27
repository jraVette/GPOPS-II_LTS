%-------------------------------------------%
% BEGIN: function particleMotionEndpoint.m %
%-------------------------------------------%
function output = fourwheelEndpoint(input)
% cost = input.auxdata.cost;
% %Objective

output.objective = -input.phase(1).finalstate(1);

end
%-------------------------------------------%
% END: function brachistochroneEndpoint.m   %
%-------------------------------------------%
