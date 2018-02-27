function daq = generateInitialDaq()
%Moving parameters to daq file
%2017-02-24 10:00 - tried to make it a sparse guess for each subsequent
%     horizon. 
%2017-02-24 10:47 - used beginning and end of last horizon + end of new
%     horizon as guess. Found that the sparse guess didn't help too much.
%     Did find that the horizon that failed. Most of the horizons
%     converged. So, going to try a updated guess if one of the meshes gave
%     a good solution
%2017-02-25 14:05 - multiple tries if horizon has one mesh w/ solution,
%     less mesh refinement
%2017-02-26 9:20 - guess from exisitng solution (took 5 mesh iterations
%     here). Then lowered mesh iters to 3.
%2017-02-26 9:52 - didn't like the file initial guess. Need to think on
%     that. Going back to 5mesh iter and arbitrary initial guess to basline
%2017-02-26 14:00 - Tried scaling again. Fixed several errors but, no
%     better than auto scaling!

%Other things to consider;
% - maybe lower mesh tolerance way down
% - guesses of subsequent horizons don't exactly match x0 conditions. Look
%   into

%% Vehicle
% vehicleDirectory = fullfile(jatecPath,'Resources/In house code/Vehicle Parameters/');
% fullVehicleFile = fullfile(vehicleDirectory,'Corvette',carFilename);
% fullVehicleFile = fullfile(vehicleDirectory,'Optimal Control Research',carFileName);
% load(fullVehicleFile);
carFilename = '2015_Corvette_C7R.mat';
% carFilename = 'LimebeerF1Car.mat';
load(carFilename);

%% Track
trackFilename = 'SebringLoopedOptFitCtrLine.mat';
refDaqFile = load('DriverA_Sebring_T7extract.mat');
% trackFilename = 'dragStrip';

track = load(trackFilename);
track = track.track;

%% Variable naming
variableNames.indepVarName = 'distance';
variableNames.stateNames   = { 'vx';'vy';'yawRate';'omegaWheel_L1';'omegaWheel_R1';'omegaWheel_L2';'omegaWheel_R2';'torqueDemand';'ey';'ePsi';'delta'};
variableNames.controlNames = {'u1';'u2'};
variableNames.units        = {'m/s';'m/s';'rad/s';'rad/s';'rad/s';'rad/s';'rad/s';'N*m';'m';'rad';'rad';'rad/s';'N*m/s';'m'};
variableNames.names        = {'v_x';'v_y';'\dot{\psi}';'\omega_{L1}';'\omega_{R1}';'\omega_{L2}';'\omega_{R2}';'T';'e_y';'e_{\psi}';'\delta';'u_1';'u_2';'s'};

%% MPC parameters
horizon                 = 250;                                             %[m] Look ahead %150m for chicane, updated based on course DOE
controlHorizon          = 10;                                              %[m] MPC update %5m for chicane, updated based on course DOE
interpolationAccuracy   = 0.25;                                            %[m] ds
horizonDecrement        = 10;                                              %[m] used to shorten horizon incase of convergence error
minimumHorizon          = 50;                                              %[m] minimum acceptable horizon
initialDistance         = refDaqFile.daq.rawData.distance.meas(1);                                            %[m] s0 %I want to start well before start/finish line                                         %Where the car gets on the track
timingDistanceStart     = initialDistance+10;                                               %Where timing starts
timingDistanceFinish    = refDaqFile.daq.rawData.distance.meas(end);
finishDistance          = timingDistanceFinish+10;
horizonRefinement       = true;

%% Scaling
scaling.length = 1;%1/(vehicle.parameter.a.meas+vehicle.parameter.b.meas);
scaling.mass   = 1;%1/(vehicle.parameter.mass.meas);
scaling.time   = 1;%sqrt(9.81/(vehicle.parameter.a.meas+vehicle.parameter.b.meas));
scaling.angle  = 1;

% scaling.length = 1/(vehicle.parameter.a.meas+vehicle.parameter.b.meas);
% scaling.mass   = 1/(vehicle.parameter.mass.meas);
% scaling.time   = sqrt(9.81/(vehicle.parameter.a.meas+vehicle.parameter.b.meas));
% scaling.angle  = 1;

scaling.velocity = scaling.length/scaling.time;
scaling.acceleration = scaling.length/scaling.time^2;
scaling.angularVelocity = scaling.angle/scaling.time;
scaling.force = scaling.mass*scaling.acceleration;
scaling.torque = scaling.force*scaling.length;
scaling.state = [scaling.velocity scaling.velocity scaling.angularVelocity scaling.angularVelocity scaling.angularVelocity scaling.angularVelocity scaling.angularVelocity scaling.torque scaling.length scaling.angle scaling.angle];
 scaling.control = [scaling.angle/scaling.time  scaling.torque/scaling.time];


%% Boundary conditions
vx0 = refDaqFile.daq.rawData.vx.meas(1);%10;
vy0 = refDaqFile.daq.rawData.vy.meas(1);%0;
r0  =refDaqFile.daq.rawData.yawRate.meas(1);
omega_front0 = vx0*(1)./vehicle.tire_front.reff.meas;
omega_rear0  = vx0*(0.0826371923756939+1)./vehicle.tire_rear.reff.meas; 
% omega_rear0  = vx0*(1)./vehicle.tire_rear.reff.meas; 
% T0 = 4110;
% T0 = 2393.5092066619;
% T0 = 0;
T0 = 2000;

ey0 = refDaqFile.daq.rawData.ey.meas(1);
ePsi0 = refDaqFile.daq.rawData.ePsi.meas(1);
delta0 = refDaqFile.daq.rawData.wheelSteerAngle_L1.meas(1)*myConstants.deg2rad;

x0 = [vx0 vy0 r0 omega_front0 omega_front0 omega_rear0 omega_rear0 T0 ey0 ePsi0 delta0];


%% Guess
%Open loop ODE guess
% sGuess = (s0:interpolationAccuracy:sf)';
% uGuess = 0*ones(size(sGuess));
% guessDaq.header.setup = setup; %Mimic how we'll have the daq in the final form for OCP
% guessDaq = generateGuessDaq(sGuess,x0,uGuess,guessDaq);
% setup.guess.phase.time    = writeDaqChannelsToMatrix(guessDaq,'selectedChannels',variableNames.indepVarName);
% setup.guess.phase.state   = writeDaqChannelsToMatrix(guessDaq,'selectedChannels',variableNames.stateNames);
% setup.guess.phase.control = writeDaqChannelsToMatrix(guessDaq,'selectedChannels',variableNames.controlNames);
% setup.guess.phase.integral     = 0;

%Existing file as guess
% guessFile = load('2018-01-24_11_17_25_solutionSnapshot_Horizon001.mat');
% guessFile = load('2018-02-06_16_02_23_solutionSnapshot_Horizon001_corvette.mat');
% guessFile = load('guess2018-02-20_17_11_56_solutionSnapshot_Horizon001.mat');
% guessFile = load('guess2018-02-26_09_06_15_solutionSnapshot_Horizon001.mat');
% setup.guess.phase.time = guessFile.segDaq.gpopsOutput.result.solution.phase.time;
% setup.guess.phase.state = [guessFile.segDaq.gpopsOutput.result.solution.phase.state];
% setup.guess.phase.control = [guessFile.segDaq.gpopsOutput.result.solution.phase.control];
% setup.guess.phase.integral = guessFile.segDaq.gpopsOutput.result.solution.phase.integral;

 
%Near arbitrary initial guess
setup.guess.phase.time    = [0; horizon];
setup.guess.phase.state   = [x0; x0];
setup.guess.phase.control = zeros(2,length(variableNames.controlNames));
setup.guess.phase.control(1,2) = 20;
setup.guess.phase.integral     = 0;
% 
setup.guess.phase.time = setup.guess.phase.time*scaling.length;
setup.guess.phase.state = bsxfun(@times,setup.guess.phase.state,scaling.state);
setup.guess.phase.control = bsxfun(@times,setup.guess.phase.control,scaling.control);




%% Bounds
vxLb        = 20;                                                          %Original bound
vxUb        = 120;%69.9240505593388;                                   %Original Bound
vyMax       = 5;                                                          %Orignal bounds
rMax        = 55*myConstants.deg2rad;                                      %Orignal bound 45 deg/s
omegaLb     = vxLb/vehicle.tire_front.reff.meas;                           %Just using the reff of the front should be sufficient
omegaUb     = vxUb/vehicle.tire_front.reff.meas;
TMax        = 5000;
TRate       = 1000;                                                        % N*m/s
eyMax       = 12.2622/2;                                                   %Road width constraint
ePsiMax     = 25*myConstants.deg2rad;                                      %Pevious solutions said this was bounded by [-25, 25]
deltaMax    = 18*myConstants.deg2rad;
deltaRate   = 100*myConstants.deg2rad;

setup.bounds.phase.initialtime.lower  = 0*scaling.length; 
setup.bounds.phase.initialtime.upper  = 0*scaling.length;
setup.bounds.phase.finaltime.lower    = horizon*scaling.length;
setup.bounds.phase.finaltime.upper    = horizon*scaling.length;
setup.bounds.phase.initialstate.lower = x0.*scaling.state;
setup.bounds.phase.initialstate.upper = x0.*scaling.state;
setup.bounds.phase.state.lower        = [vxLb  -vyMax -rMax omegaLb omegaLb omegaLb omegaLb -TMax -eyMax -ePsiMax -deltaMax].*scaling.state; 
setup.bounds.phase.state.upper        = [vxUb   vyMax  rMax omegaUb omegaUb omegaUb omegaUb  TMax  eyMax  ePsiMax  deltaMax].*scaling.state;
setup.bounds.phase.finalstate.lower   = setup.bounds.phase.state.lower;
setup.bounds.phase.finalstate.upper   = setup.bounds.phase.state.upper;
setup.bounds.phase.control.lower      = [ -deltaRate -TRate].*scaling.control;
setup.bounds.phase.control.upper      = [  deltaRate  TRate].*scaling.control;
setup.bounds.phase.path.lower         = [0*ones(1,4) ];%, -100];
setup.bounds.phase.path.upper         = [0.8*ones(1,4) ];%,  100];
setup.bounds.phase.integral.lower     =  0;
setup.bounds.phase.integral.upper     =  1e9;

%% Auxdata - put this here so I can use bounds
setup.auxdata.variableNames             = variableNames;
setup.auxdata.scaling                   = scaling;
setup.auxdata.vehicle                   = vehicle;
setup.auxdata.track                     = track;
setup.auxdata.minTimeCost               = 1;
setup.auxdata.controlWeight             = 1e-6;%1e-3;
setup.auxdata.torqueAllocationSlope     = 10; %used in the sin(atan(.)) to seperate plus and minus. 1 seemed to work fine
setup.auxdata.currentDistance           = initialDistance;
% vxTarget = 100;   
% setup.auxdata.stageCost.targetState     = [vxTarget 0 0 omegaUb omegaUb omegaUb omegaUb 0 0 0 0];
% setup.auxdata.stageCost.scaling         = [1        1 1 0       0       0       0       0 1 1 0];
% setup.auxdata.stageCost.weight          = [1];

% setup.auxdata.terminalCost.targetState  = [vxTarget 0 0 omegaUb omegaUb omegaUb omegaUb 0 0 0];
% setup.auxdata.terminalCost.scaling      = [vxTarget vyMax rMax omegaUb omegaUb omegaUb omegaUb TMax eyMax ePsiMax];
% setup.auxdata.terminalCost.weight       = 0;


%% GPOPS Setup
setup.name                        = 'quadCar';
setup.nlp.solver                  = 'ipopt';
setup.derivatives.supplier        = 'adigator';%'adigator';%'sparseFD'; %'adigator';
setup.derivatives.derivativelevel = 'second';
setup.scales.method               = 'automatic-hybridUpdate';%'automatic-guessUpdate';'automatic-hybridUpdate';'none'
% 'automatic-bounds'       scales the problem from the user-supplied bounds on the variables
% 'automatic-guess'        scales the problem once using the initial guess of the solution supplied by the user 
% 'automatic-guessUpdate?  scales the problem from the initial guess on the first mesh and from the solution obtained on every subsequent mesh during the mesh refinement
% 'automatic-hybrid'       scales the problem from the user supplied bounds on the variables on the first mesh and from the solution obtained on the initial mesh for every subsequent mesh in the mesh refinement
% 'automatic-hybridUpdate' scales the problem from the bounds on the initial mesh and from the solution obtained on every subsequent mesh during the mesh refinement


setup.method                      = 'RPM-Integration';
setup.displaylevel                = 2;
setup.nlp.ipoptoptions.maxiterations = 1000;

setup.mesh.method       = 'hp-PattersonRao';
setup.mesh.tolerance    = 1e-2;
setup.mesh.maxiterations = 5;
nFrac = 4;
setup.mesh.phase.fraction = 1/nFrac*ones(1,nFrac);
setup.mesh.phase.colpoints = 4*ones(1,nFrac);
acceptableNlpOutpus = [0 1 ] ; 


%% DAQ File
filename = sprintf('%s_GPOPS_ShortSegStrightLine-tOpt',datestr(now,'yyyy-mm-dd_HH_MM_SS'));
shortFilename = 'tOpt';
daq.header = saveVariablesAssignedToPointInStructure('exclude',{'varargin';'vehicle';'track'},'clearVariableAfterPackage',true);
daq.vehicle = vehicle;
daq.track = track;
% daq.header.iterNumb = 1;
daq.header.path = pwd;
daq.header = addNotesToDaqFile(daq.header,sprintf('File setup %s to setup for MPC LTS',datestr(now,'yyyy-mm-dd_HH_MM_SS')));



