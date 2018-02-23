function guessDaq = generateGuessDaq(s,x0,uGuess,guessDaq)
%This function will run open loop dynamics of the ODE around an initial
%guess.
%INPTUS:
%     s - discritized vector of the independant variable space size 1xt
%     x0 - initial states nx1
%     uGuess - vector size of mxt w/ the control
%     vehicle - standard vehicle defnition structure
%     track - standard track structure
%     guessDaq, daq with a .auxdata field
%OUTPUS:
%     guessDaq - open loop daq simulation
%
%Creation: 20 Dec 2017 - Jeff Anderson
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

if nargin == 0
    guessDaq = generateInitialDaq;
    inputDaq = load('guess2018-02-20_17_11_56_solutionSnapshot_Horizon001.mat');
    getChannelDataFromDaqFile(inputDaq.segDaq,{'distance';'u1';'u2'});
    s = distance.*guessDaq.header.scaling.length;
    uGuess = bsxfun(@times,[u1 u2],guessDaq.header.scaling.control);

    x0 = guessDaq.header.x0;
end

setup = guessDaq.header.setup; %Save it out to make calling shorter and we'll have to re-insert it in the parsed solution

sSpan = [s(1) s(end)];

options = odeset('MaxStep',0.1);
[sSol,xSol] = ode45(@(sInput,x)fourwheelMatlabODE(sInput,x,uGuess,s,setup.auxdata), sSpan, x0,options); 

     
%Resample to requested input
SSOL = s;
for i = 1:size(xSol,2)
    XSOL(:,i) = interp1(sSol,xSol(:,i),SSOL);
end
SSOL = SSOL./guessDaq.header.scaling.length;
XSOL = bsxfun(@times,XSOL,1./guessDaq.header.scaling.state);
uGuess = bsxfun(@times,uGuess,1./guessDaq.header.scaling.control);

guessDaq = parseMatlabOdeOutput(SSOL,XSOL,uGuess,...
    setup.auxdata.variableNames.indepVarName,...
    setup.auxdata.variableNames.stateNames,...
    setup.auxdata.variableNames.controlNames,...
    setup.auxdata.variableNames.units,...
    setup.auxdata.variableNames.names);
guessDaq.header.setup = setup; %Put it back in

%Mimic gpops output for comparison
solutionField = 'interpsolution'; %field in the gpops results output.results.(solutionField) originally just solution but looks like the interpsol is more what I need.
gpopsOutput.result.(solutionField).phase.time = SSOL;
gpopsOutput.result.(solutionField).phase.state = XSOL;
gpopsOutput.result.(solutionField).phase.control = uGuess;
gpopsOutput = addNotesToDaqFile(gpopsOutput,'Just added this to mimic gpops output for calculating algebraic states');
guessDaq.gpopsOutput = gpopsOutput; 
guessDaq = calculateAlgebraicStates(guessDaq);

