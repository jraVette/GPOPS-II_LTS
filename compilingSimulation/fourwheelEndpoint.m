%-------------------------------------------%
% BEGIN: function particleMotionEndpoint.m %
%-------------------------------------------%
function output = fourwheelEndpoint(input)
costTime = input.auxdata.costTime;
costVx = input.auxdata.costVx;

% %Objective
integral = input.phase(1).integral;
% terminalCosts = sum(((input.auxdata.terminalCost.targetState - input.phase.finalstate)./input.auxdata.terminalCost.scaling).^2)/length(input.phase.finalstate) ;
% output.objective = integral + terminalCosts*input.auxdata.terminalCost.weight;
vxf = input.phase.finalstate(1);
tf = input.phase.finalstate(end);
% J_tf =    tf + input.auxdata.regularizationCost*integral(1);
J_vx =  -vxf + input.auxdata.regularizationCostVx*integral;
J_tf = tf + input.auxdata.regularizationCost*integral; 

% output.objective = cost*J_vx/input.auxdata.vxCostScale + (1-cost)*J_tf/input.auxdata.tCostScale;

% output.objective = 0*J_tf + 1*J_vx;
% output.objective = 1*J_vx + 0*123;
output.objective = costVx*J_vx + costTime*J_tf;
% output.objective =  costTime*J_tf;
% output.eventgroup.event = [input.phase.finalstate - input.phase.initialstate];

end
%-------------------------------------------%
% END: function brachistochroneEndpoint.m   %
%-------------------------------------------%
