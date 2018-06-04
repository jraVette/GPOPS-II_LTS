function daq = generateInitialDaq(varargin)
%OCP over fullfile

defaults = {'loadGuess',true};
setDefaultsForVarargin(defaults,varargin)

%% Vehicle
% vehicleDirectory = fullfile(jatecPath,'Resources/In house code/Vehicle Parameters/');
% fullVehicleFile = fullfile(vehicleDirectory,'Corvette',carFilename);
% fullVehicleFile = fullfile(vehicleDirectory,'Optimal Control Research',carFileName);
% load(fullVehicleFile);
 carFilename = '2015_Corvette_C7R.mat'
% carFilename = 'LimebeerF1Car.mat';
load(carFilename);

refDaqFile = 'driverB_EstInputs_T17.mat'
refDaq = load(refDaqFile);
refDaq = refDaq.daq;

%% Track
%trackFilename = 'PathInfoChicaneStraightsBeforeAndAfter.mat';
%trackFilename = 'dragStrip.mat';
 trackFilename = 'SebringLoopedOptFitCtrLine.mat'
track = load(trackFilename);
track = track.track;

%% Variable naming
variableNames.indepVarName = 'distance';
variableNames.stateNames   = { 'vx';'vy';'yawRate';'ey';'ePsi';'time'};
variableNames.controlNames = {'delta';'kappa_L1';'kappa_R1';'kappa_L2';'kappa_R2';'fz_L1';'fz_R1';'fz_L2';'fz_R2'};
variableNames.units        = {'m/s';'m/s';'rad/s';     'm';' rad';'s';      'rad';  '';        '';        '';        '';        'N';    'N';    'N';     'N'; 'm'};
variableNames.names        = {'v_x';'v_y';'\dot{\psi}';'e_y';'e_{\psi}';'time';'\delta';'\kappa_{L1}';'\kappa_{R1}';'\kappa_{L2}';'\kappa_{R2}';'fz_{L1}';'fz_{R1}';'fz_{L2}';'fz_{R2}';'distance'};


clear s c %don't need them

%% MPC parameters
horizon                 = 300;                                             %[m] Look ahead %150m for chicane, updated based on course DOE
controlHorizon          = 60;                                              %[m] MPC update %5m for chicane, updated based on course DOE
interpolationAccuracy   = 0.25;                                            %[m] ds
horizonDecrement        = 5;                                              %[m] used to shorten horizon incase of convergence error
minimumHorizon          = 50;                                              %[m] minimum acceptable horizon
initialDistance         = 4565;                                            %[m] s0 %I want to start well before start/finish line                                         %Where the car gets on the track
timingDistanceStart     = 4565;                                               %Where timing starts
timingDistanceFinish    = 5859;
finishDistance          = timingDistanceFinish+10;
horizonRefinement       = false;


%% Setup switching
s                       = [-100 1*10000]';                                 %just big numbers so it spans the track
c                       = [0 0]';                                          %time optimal switching
switchingDaq.rawData.distance = createDaqChannelData(s,'m','Distance');
switchingDaq.rawData.switching = createDaqChannelData(c,'','Switching');

s = initialDistance:controlHorizon:finishDistance;
     %'time cost vx cost'},...
c0  = [ 1        0];
lb0 = [ 0        0 ];
ub0 = [ 1        1 ];
c  = repmat(c0,length(s),[]);
lb = repmat(lb0,length(s),[]);
ub = repmat(ub0,length(s),[]);

switchingDaq.rawData.distance = createDaqChannelData(s,'m','Distance');
switchingDaq.rawData.switching = createDaqChannelData(c,'','Switching');
switchingDaq.rawData.lb = createDaqChannelData(lb,'','Switching Lower Bounds');
switchingDaq.rawData.ub = createDaqChannelData(ub,'','Switching Upper Bounds');
clear s c c0 lb ub lb0 ub0%don't need them

%% Scaling
scaling.length = 1;%1/(vehicle.parameter.a.meas+vehicle.parameter.b.meas);
scaling.mass   = 1;%1/(vehicle.parameter.mass.meas);
scaling.time   = 1;%sqrt(scaling.length*9.81);
scaling.angle  = 1;

scaling.velocity = scaling.length/scaling.time;
scaling.acceleration = scaling.length/scaling.time^2;
scaling.angularVelocity = scaling.angle/scaling.time;
scaling.force = scaling.mass*scaling.acceleration;
scaling.torque = scaling.force*scaling.length;
scaling.state = [scaling.velocity scaling.velocity scaling.angularVelocity scaling.length scaling.angle scaling.time];
scaling.control = [scaling.angle 1 1 1 1 scaling.force scaling.force scaling.force scaling.force];

%% Boundary conditions
ind  = findNearestPoint(initialDistance,refDaq.rawData.distance.meas);
x0Daq = assembleNewDaqAtIndicies(ind,refDaq);
x0   = writeDaqChannelsToMatrix(x0Daq,'selectedChannels',variableNames.stateNames);
clear x0Daq

x0(end) = 0;
x0 = x0.*scaling.state;


%% Bounds
vxLb        = 5;                                                           %Original bound
vxUb        = 160;%69.9240505593388;                                       %Original Bound
vyMax       = 10;                                                          %Orignal bounds
rMax        = 70*myConstants.deg2rad;                                      %Orignal bound 45 deg/s
eyMax       = 6.1311;                                                           %Road width constraint
ePsiMax     = 25*myConstants.deg2rad;                                      %Pevious solutions said this was bounded by [-25, 25]
deltaMax    = 30*myConstants.deg2rad;
tMax        = 30;
kappaMin = -0.09;
kappaMax = 0.09;
fzMin = -30000;
fzMax = 0;

setup.bounds.phase.initialtime.lower  = 0*scaling.length; 
setup.bounds.phase.initialtime.upper  = 0*scaling.length;
setup.bounds.phase.finaltime.lower    = horizon*scaling.length;
setup.bounds.phase.finaltime.upper    = horizon*scaling.length;
setup.bounds.phase.initialstate.lower = x0.*scaling.state;
setup.bounds.phase.initialstate.upper = x0.*scaling.state;
setup.bounds.phase.state.lower        = [vxLb  -vyMax -rMax -eyMax -ePsiMax 0   ].*scaling.state; 
setup.bounds.phase.state.upper        = [vxUb   vyMax  rMax  eyMax  ePsiMax tMax].*scaling.state;
setup.bounds.phase.finalstate.lower   = setup.bounds.phase.state.lower.*scaling.state;
setup.bounds.phase.finalstate.upper   = setup.bounds.phase.state.upper.*scaling.state;
setup.bounds.phase.control.lower      = [ -deltaMax 0        0        kappaMin kappaMin fzMin fzMin fzMin fzMin].*scaling.control;
setup.bounds.phase.control.upper      = [  deltaMax kappaMax kappaMax kappaMax kappaMax fzMax fzMax fzMax fzMax].*scaling.control;
setup.bounds.phase.path.lower         = [0*ones(1,6) 0 ];
setup.bounds.phase.path.upper         = [0*ones(1,6) 1];
setup.bounds.phase.integral.lower     =  [0  ];
setup.bounds.phase.integral.upper     =  [1e8];

%% Auxdata - put this here so I can use bounds
setup.auxdata.variableNames             = variableNames;
setup.auxdata.scaling                   = scaling;
setup.auxdata.vehicle                   = vehicle;
setup.auxdata.track                     = track;
setup.auxdata.minTimeCost               = 1;
setup.auxdata.controlWeight             = 1e-6;%1e-3;
setup.auxdata.torqueAllocationSlope     = 10; %used in the sin(atan(.)) to seperate plus and minus. 1 seemed to work fine
setup.auxdata.currentDistance           = initialDistance;
setup.auxdata.controlCost               = [1.5e1 2e-1 2e-1 2e-1 2e-1 2e-11 2e-11 2e-11 2e-11];
setup.auxdata.regularizationCost        = 1e-0;%was 1e-2

setup.auxdata.controlCostVx             = [1e2 2e0 2e0 2e0 2e0 2e-11 2e-11 2e-11 2e-11];
setup.auxdata.regularizationCostVx      = 1e-0;%was 1e-2
setup.auxdata.tCostScale                = 1;%30;
setup.auxdata.vxCostScale               = 1;%-70;
setup.auxdata.engMult                   = 0.9;
setup.auxdata.muMultX                   = 0.58;
setup.auxdata.muMultY                   = 1.6;
setup.auxdata.costTime                  = 1; %Note gpopsMPC will overwrite this, need it to compile adigator
setup.auxdata.costVx                    = 1; %Note gpopsMPC will overwrite this, need it to compile adigator



%% GPOPS Setup
setup.name                        = 'quadCar';
setup.nlp.solver                  = 'ipopt';
setup.derivatives.supplier        = 'adigator';%'adigator';%'adigator';%'sparseFD'; %'adigator';
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
setup.mesh.tolerance    = 1e-3;
setup.mesh.maxiterations = 4;
nFrac = 4;
setup.mesh.phase.fraction = 1/nFrac*ones(1,nFrac);
setup.mesh.phase.colpoints = 4*ones(1,nFrac);
setup.mesh.colpointsmax = 6;
setup.mesh.colpointsmin = 2;
setup.mesh
acceptableNlpOutputs = [0 1 ] ; 
nextGuessType = 'exactPreviousGuess'; % basedOffPrevious  basedOffPreviousStartAndEndPoints  none %PRETTY SURE UNUSED DEPRECATE



%% Guess
if loadGuess
    g = 9.81;
    m = vehicle.parameter.mass.meas;
    a = vehicle.parameter.a.meas;
    b = vehicle.parameter.b.meas;
    wf = -b/(a+b)*m*g/2;
    wr = -a/(a+b)*m*g/2;

    %Near arbitrary initial guess
    setup.guess.phase.time    = [0; finishDistance-initialDistance];
    setup.guess.phase.state   = [x0; [x0(1:end-1) 5]];
    setup.guess.phase.control = ...
        [0 0.00 0.00 0.00 0.00 wf wf wr wr
         0 0.00 0.00 0.00 0.00 wf wf wr wr];
    setup.guess.phase.integral     = [0];
    
    setup.guess.phase.time = setup.guess.phase.time*scaling.length;
    setup.guess.phase.state = bsxfun(@times,setup.guess.phase.state,scaling.state);
    setup.guess.phase.control = bsxfun(@times,setup.guess.phase.control,scaling.control);
end

%% Ga options
populationSize = 96; 
nIteratesToSavePerGeneration = 10; %Leave empty for all
timeToGiveUpOnSim = 100*60*60;
nonConvergentCost = 1000;
simFinished = false; %Flag for stats
conv = nan; %Flag for stats
saveDiaryFiles = false;
cleanUpRunningDirecotry = true;

gaCost.independentVarChannel = 'distance';
    %Ch name    %Weight Q   %Scaling S
tempCost = {  
    'delta'     0           1.142035744625003e-01
    'ePsi'      0           2.590834982889016e+00
    'ey'        2           7.369512523277982e+03
    'fz_L1'     0           2.507200546069650e+08
    'fz_L2'     0           2.978229664338260e+08
    'fz_R1'     0           7.521499820961577e+07
    'fz_R2'     0           2.771859743041586e+08
    'kappa_L1'  0            3.237185060229598e-02
    'kappa_L2'  0            1.126754617810116e-01
    'kappa_R1'  0            1.550567448709579e-01
    'kappa_R2'  0            1.885937567533368e-01
    'time'      0           1.189521250818846e+07
    'vx'        1           1.494813341729434e+04
    'vy'        0           2.444859968954880e+02
    'yawRate'   0           1.760359618786300e+00
};    
gaCost.channelsToCompare = tempCost(:,1);
gaCost.weighting = cell2matrix(tempCost(:,2));
gaCost.scaling = cell2matrix(tempCost(:,3));
clear tempCost



%% DAQ File
% filename = sprintf('%s_GPOPS_ShortSegStrightLine-tOpt',datestr(now,'yyyy-mm-dd_HH_MM_SS'));
filename = 'daqFile';
shortFilename = 'gaIterate';
daq.header = saveVariablesAssignedToPointInStructure('exclude',{'varargin';'vehicle';'track'; 'refDaq'},'clearVariableAfterPackage',true);
daq.vehicle = vehicle;
daq.track = track;
% daq.header.iterNumb = 1;
daq.header.path = pwd;
daq.header = addNotesToDaqFile(daq.header,sprintf('File setup %s to setup for MPC LTS',datestr(now,'yyyy-mm-dd_HH_MM_SS')));



