function daq = generateInitialDaq()
%Moving parameters to daq file

%% Vehicle
vehicleDirectory = fullfile(jatecPath,'Resources/In house code/Vehicle Parameters/');
% carFilename = '2015_Corvette_C7R.mat';
% fullVehicleFile = fullfile(vehicleDirectory,'Corvette',carFilename);
carFileName = 'LimebeerF1Car.mat';
fullVehicleFile = fullfile(vehicleDirectory,'Optimal Control Research',carFileName);
load(fullVehicleFile);

%% Track
trackFilename = 'PathInfoChicaneStraightsBeforeAndAfter.mat';
trackFilename = 'dragStrip.mat';

track = load(trackFilename);
track = track.track;

%% Variable naming
variableNames.indepVarName = 'distance';
variableNames.stateNames   = { 'vx';'vy';'yawRate';'omegaWheel_L1';'omegaWheel_R1';'omegaWheel_L2';'omegaWheel_R2';'torqueDemand';'ey';'ePsi'};
variableNames.controlNames = {'u2'};
variableNames.units        = {'m';'m/s';'m/s';'rad/s';'rad/s';'rad/s';'rad/s';'rad/s';'N*m';'m';'rad';'N*m/s'};
variableNames.names        = {'Time';'Vx';'Vy';'Yaw Rate';'Wheel Speed Left Front';'Wheel Speed Right Front';'Wheel Speed Left Rear';'Wheel Speed Right Rear';'Torque Demand';'Lateral Deviation';'Heading Deviation';'Torque Demand Rate'};

%% MPC parameters
horizon                 = 200;                                             %[m] Look ahead %150m for chicane, updated based on course DOE
controlHorizon          = 10;                                              %[m] MPC update %5m for chicane, updated based on course DOE
interpolationAccuracy   = 0.25;                                            %[m] ds
horizonDecrement        = 10;                                              %[m] used to shorten horizon incase of convergence error
minimumHorizon          = 50;                                              %[m] minimum acceptable horizon
initialDistance         = -200;                                            %[m] s0 %I want to start well before start/finish line                                         %Where the car gets on the track
timingDistanceStart     = 0;                                               %Where timing starts
timingDistanceFinish    = track.finishDistance;
finishDistance          = timingDistanceFinish+10;
horizonRefinement       = true;


%% Boundary conditions
s0 = initialDistance;
sf = initialDistance+horizon;
vx0 = 10;
vy0 = 0;
r0  = 0;
omega_front0 = vx0*(1)./vehicle.tire_front.reff.meas;
omega_rear0  = vx0*(0.0826371923756939+1)./vehicle.tire_rear.reff.meas; 
% omega_rear0  = vx0*(1)./vehicle.tire_rear.reff.meas; 
% T0 = 4110;
T0 = 2393.5092066619;
% T0 = 0;

ey0 = 0;
ePsi0 = 0;

x0 = [vx0 vy0 r0 omega_front0 omega_front0 omega_rear0 omega_rear0 T0 ey0 ePsi0];



%% Guess
% sGuess = (s0:interpolationAccuracy:sf)';
% uGuess = 0*ones(size(sGuess));
% guessDaq.header.setup = setup; %Mimic how we'll have the daq in the final form for OCP
% guessDaq = generateGuessDaq(sGuess,x0,uGuess,guessDaq);
% setup.guess.phase.time    = writeDaqChannelsToMatrix(guessDaq,'selectedChannels',variableNames.indepVarName);
% setup.guess.phase.state   = writeDaqChannelsToMatrix(guessDaq,'selectedChannels',variableNames.stateNames);
% setup.guess.phase.control = writeDaqChannelsToMatrix(guessDaq,'selectedChannels',variableNames.controlNames);
% setup.guess.phase.integral     = 0;

%Load file
guessFile = load('2018-01-24_11_17_25_solutionSnapshot_Horizon001.mat');
setup.guess.phase.time = guessFile.segDaq.gpopsOutput.result.solution.phase.time;
setup.guess.phase.state = guessFile.segDaq.gpopsOutput.result.solution.phase.state;
setup.guess.phase.control = guessFile.segDaq.gpopsOutput.result.solution.phase.control;
setup.guess.phase.integral = guessFile.segDaq.gpopsOutput.result.solution.phase.integral;

%Near arbitrary initial guess
% setup.guess.phase.time    = [s0; sf];
% setup.guess.phase.state   = [x0; x0];
% setup.guess.phase.control = zeros(2,length(variableNames.controlNames));
% setup.guess.phase.control(1,1) = 0.5;
% setup.guess.phase.integral     = 0;




%% Bounds
vxLb        = 0;                                                           %Original bound
vxUb        = 150;%69.9240505593388;                                       %Original Bound
vyMax       = 10;                                                        %Orignal bounds
rMax        = 55*myConstants.deg2rad;                                    %Orignal bound 45 deg/s
omegaLb     = vxLb/vehicle.tire_front.reff.meas;                           %Just using the reff of the front should be sufficient
omegaUb     = vxUb/vehicle.tire_front.reff.meas;
TMax        = 5000;
TRate       = 100*1000/5000;                                                 % N*m/s
eyMax       = 5;                                                       %Road width constraint
ePsiMax     = 25*myConstants.deg2rad;                                  %Pevious solutions said this was bounded by [-25, 25]

setup.bounds.phase.initialtime.lower  = s0; 
setup.bounds.phase.initialtime.upper  = s0;
setup.bounds.phase.finaltime.lower    = sf; 
setup.bounds.phase.finaltime.upper    = sf;
setup.bounds.phase.initialstate.lower = x0;
setup.bounds.phase.initialstate.upper = x0;
setup.bounds.phase.state.lower        = [vxLb  -vyMax -rMax omegaLb omegaLb omegaLb omegaLb -TMax -eyMax -ePsiMax]; 
setup.bounds.phase.state.upper        = [vxUb   vyMax  rMax omegaUb omegaUb omegaUb omegaUb  TMax  eyMax  ePsiMax];
setup.bounds.phase.finalstate.lower   = setup.bounds.phase.state.lower;
setup.bounds.phase.finalstate.upper   = setup.bounds.phase.state.upper;
setup.bounds.phase.control.lower      = [ -1];
setup.bounds.phase.control.upper      = [  1];
setup.bounds.phase.path.lower         = [-0.2*ones(1,4)];
setup.bounds.phase.path.upper         = [ 0.2*ones(1,4)];
setup.bounds.phase.integral.lower     =  0;
setup.bounds.phase.integral.upper     =  1e9;

%% Auxdata - put this here so I can use bounds
setup.auxdata.variableNames             = variableNames;
setup.auxdata.vehicle                   = vehicle;
setup.auxdata.track                     = track;
setup.auxdata.minTimeCost               = 1;
setup.auxdata.controlWeight             = 1e-6;%1e-3;
setup.auxdata.torqueAllocationSlope     = 10; %used in the sin(atan(.)) to seperate plus and minus. 1 seemed to work fine
vxTarget = 100;
setup.auxdata.stageCost.targetState     = [vxTarget 0 0 omegaUb omegaUb omegaUb omegaUb 0 0 0];
setup.auxdata.stageCost.scaling         = [1        1 1 0       0       0       0       0 1 1];
setup.auxdata.stageCost.weight          = [1];

setup.auxdata.terminalCost.targetState  = [vxTarget 0 0 omegaUb omegaUb omegaUb omegaUb 0 0 0];
setup.auxdata.terminalCost.scaling      = [vxTarget vyMax rMax omegaUb omegaUb omegaUb omegaUb TMax eyMax ePsiMax];
setup.auxdata.terminalCost.weight       = 0;


%% GPOPS Setup
setup.name                        = 'quadCar';
setup.nlp.solver                  = 'ipopt';
setup.derivatives.supplier        = 'adigator';%'adigator';%'sparseFD'; %'adigator';
setup.derivatives.derivativelevel = 'second';
setup.scales.method               = 'automatic-hybridUpdate';
setup.method                      = 'RPM-Integration';
setup.displaylevel                = 2;
setup.nlp.ipoptoptions.maxiterations = 1000;

setup.mesh.method       = 'hp-PattersonRao';
setup.mesh.tolerance    = 1e-3;
setup.mesh.maxiterations = 10;
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


%% Create adigator files
setup = daq.header.setup;
setup.functions.continuous        = @fourwheelContinuous;
setup.functions.endpoint          = @fourwheelEndpoint;
adigatorGenFiles4gpops2(setup);

