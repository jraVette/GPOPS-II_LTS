function guessDaq = generateGuessDaq(s,x0,uGuess,vehicle,track,variableNames)
%This function will run open loop dynamics of the ODE around an initial
%guess.
%INPTUS:
%     s - discritized vector of the independant variable space size 1xt
%     x0 - initial states nx1
%     uGuess - vector size of mxt w/ the control
%     vehicle - standard vehicle defnition structure
%     track - standard track structure
%     varaible names - variable name structure
%OUTPUS:
%     guessDaq - open loop daq simulation
%
%Creation: 20 Dec 2017 - Jeff Anderson
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


auxdata.controlWeight = 1e-3;
auxdata.vehicle       = vehicle;
auxdata.track         = track;
auxdata.indepVarName  = variableNames.indepVarName;
auxdata.stateNames    = variableNames.stateNames;
auxdata.controlNames  = variableNames.controlNames;
auxdata.units         = variableNames.units;
auxdata.names         = variableNames.names;

sSpan = [s(1) s(end)];

options = odeset('MaxStep',0.1);
[sSol,xSol] = ode45(@(sInput,x)fourwheelMatlabODE(sInput,x,uGuess,s,auxdata), sSpan, x0,options);

     
%Resample to requested input
SSOL = s;
for i = 1:size(xSol,2)
    XSOL(:,i) = interp1(sSol,xSol(:,i),SSOL);
end


guessDaq = parseMatlabOdeOutput(SSOL,XSOL,uGuess,auxdata.indepVarName,...
    auxdata.stateNames,...
    auxdata.controlNames,...
    auxdata.units,...
    auxdata.names);
guessDaq.auxdata = auxdata;


%Mimic gpops output for comparison
solutionField = 'interpsolution'; %field in the gpops results output.results.(solutionField) originally just solution but looks like the interpsol is more what I need.
gpopsOutput.result.(solutionField).phase.time = SSOL;
gpopsOutput.result.(solutionField).phase.state = XSOL;
gpopsOutput.result.(solutionField).phase.control = uGuess;
gpopsOutput = addNotesToDaqFile(gpopsOutput,'Just added this to mimic gpops output for calculating algebraic states');
guessDaq.gpopsOutput = gpopsOutput; 
guessDaq.vehicle = vehicle;
guessDaq = calculateAlgebraicStates(guessDaq);

