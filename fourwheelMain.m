
%Load Vehicle 
vehicleDirectory = fullfile(jatecPath,'Resources/In house code/Vehicle Parameters/');
carFilename = '2015_Corvette_C7R.mat';
fullVehicleFile = fullfile(vehicleDirectory,'Corvette',carFilename);
% carFileName = 'LimebeerF1Car.mat';
% fullVehicleFile = fullfile(vehicleDirectory,'Optimal Control Research',carFileName);
load(fullVehicleFile);

%Names
indepVarName = 'time';
stateNames = { 'vx';'omegaWheel_L1';'omegaWheel_R1';'omegaWheel_L2';'omegaWheel_R2';'torqueDemand'};
controlNames = {'u2'};
units =      {'s';'m/s';'rad/s';'rad/s';'rad/s';'rad/s';'N*m';'N*m/s'};
names = {'Time';'Vx';'Wheel Speed Left Front';'Wheel Speed Right Front';'Wheel Speed Left Rear';'Wheel Speed Right Rear';'Torque Demand';'Torque Demand Rate'};

% indepVarName = 'time';
% stateNames = { 'vx';'omegaWheel_L1'};%,'torqueDemand'};
% controlNames = {'T'};
% units =      {'s';'m/s';'rad/s';'N*m'};
% names = {'Time';'Vx';'Wheel Speed Left Front';'Torque'};

% stateNames = { 'vx'};%,'torqueDemand'};
% controlNames = {'slipRatio'};
% units =      {'s';'m/s';''};
% names = {'Time';'Vx';'Slip Ratio'};


%Initalizaiton
s0          = 0;                                             %[m] s
sf          = 10;

%% Guess
%U Guess
    s0          = 0;                                             %[m] s
    sf          = 10;
    s = s0:1:sf;
    u = 0*ones(size(s));
    vx0 = 10;
    omega_front0 = vx0*(1)./vehicle.tire_front.reff.meas;
    omega_rear0 = vx0*(0.0826371924922622+1)./vehicle.tire_rear.reff.meas;
%     T0 = 0.134334825483566*5000;
%     T0 = 2015;
    T0 = 4097.06566624625;
    x0 = [vx0 omega_front0 omega_front0 omega_rear0 omega_rear0 T0];
%     x0 = vx0;
    auxdata.vehicle = vehicle;
    auxdata.indepVarName = indepVarName;
    auxdata.stateNames = stateNames;
    auxdata.controlNames = controlNames;
    auxdata.units = units;
    auxdata.names = names;
    
    guessDaq = runDoubleTrackMatlabOde(s,x0,u,auxdata);
    guessDaq = calculateAlgebraicStates(guessDaq);
    close all
plotDaqChannelsAtEachWheelPosition('time','slipRatio',guessDaq)
% guessDaq = load('snapshot.mat');
% guessDaq = guessDaq.daq;

getChannelDataFromDaqFile(guessDaq,{'s', 'time'; 'u', 'u2'})    
stateGuess = writeDaqChannelsToMatrix(guessDaq,'selectedChannels',{'vx','omegaWheel_L1','omegaWheel_R1','omegaWheel_L2','omegaWheel_R2','torqueDemand'});
% stateGuess = writeDaqChannelsToMatrix(guessDaq,'selectedChannels',{'vx','omegaWheel_L1'});%,'torqueDemand'});
% stateGuess = writeDaqChannelsToMatrix(guessDaq,'selectedChannels',{'vx'});%,'torqueDemand'});

timeGuess               = [rowVector(s)];   
x0                      = stateGuess(1,:);
controlGuess            = [rowVector(u)];




%Deal with bounds
vxLb        = 0;                                                       %Original bound
vxUb        = 150;%69.9240505593388;                                        %Original Bound
% vyMax       = 10;                                                      %Orignal bounds
% rMax        = 55*myConstants.deg2rad;                                  %Orignal bound 45 deg/s
% tLb         = 0;
% tUb         = 15;%4.609395789295020;                                       %Oringal bounds
omegaLb     = vxLb/vehicle.tire_front.reff.meas;
omegaUb     = vxUb/vehicle.tire_front.reff.meas;
% omegaLb = -1e5;
% omegaUb = 1e5;
TMax        = 5000;
TRate       = 100*1000/5000;                                                 % N*m/s

bounds.lbX              = [vxLb    omegaLb omegaLb omegaLb omegaLb   -TMax]; 
bounds.ubX              = [vxUb    omegaUb  omegaUb omegaUb omegaUb  TMax];
% bounds.lbX              = [vxLb    omegaLb];%  -TMax]; 
% bounds.ubX              = [vxUb    omegaUb];%   TMax];
% bounds.lbX              = [vxLb  ];%  -TMax]; 
% bounds.ubX              = [vxUb  ];%   TMax];

% bounds.lbU              = [ -TMax];
% bounds.ubU              = [  TMax];
bounds.lbU              = [ -1];
bounds.ubU              = [ 1];

bounds.pathLower        = [-0.2*ones(1,4)];
bounds.pathUpper        = [ 0.2*ones(1,4)];
auxdata.bounds = bounds;
% bounds.phase.integral.lower     = -1e6;
% bounds.phase.integral.upper     = 1e6;
% clear ePsiMax eyMax vxLb vxUb myMax rMax tLb tUb omegaLb omegaUb  TMax  TRate DAQ


gpopsOptions.acceptableNlpOutpus = [0 1 ] ; 
gpopsOptions.repeatNonConvergentMpcSolNTimes = 1;
gpopsOptions.specifyMeshIterationSolution = 'auto'; %Choose 'auto' for the default or an integer for the mesh number

gpopsOptions.mesh.method       = 'hp-PattersonRao';
gpopsOptions.mesh.tolerance    = 1e-5;
gpopsOptions.mesh.maxiterations = 5;
nFrac = 10;
gpopsOptions.mesh.phase.fraction = 1/nFrac*ones(1,nFrac);
gpopsOptions.mesh.phase.colpoints = 3*ones(1,nFrac);
% gpopsOptions.mesh.colpointsmin = 7;
% gpopsOptions.mesh.colpointsmax = 20;

gpopsOptions.setup.name                        = 'quadCar';
% gpopsOptions.setup.functions.continuous        = @fourwheelContinuous;
% gpopsOptions.setup.functions.endpoint          = @fourwheelEndpoint;
gpopsOptions.setup.nlp.solver                  = 'ipopt';
gpopsOptions.setup.derivatives.supplier        = 'adigator';%'adigator';%'sparseFD'; %'adigator';%
gpopsOptions.setup.derivatives.derivativelevel = 'second';
gpopsOptions.setup.scales.method               = 'automatic-hybridUpdate';
gpopsOptions.setup.method                      = 'RPM-Differentiation';
gpopsOptions.setup.displaylevel                = 2;





%Set up daq
filename = sprintf('%s_GPOPS_ShortSegStrightLine-tOpt',datestr(now,'yyyy-mm-dd_HH_MM_SS'));
shortFilename = 'tOpt';
daq.header = saveVariablesAssignedToPointInStructure('exclude',{'varargin';'vehicle'},'clearVariableAfterPackage',false);
daq.vehicle = vehicle;
daq.header.iterNumb = 1;
daq.header.path = pwd;
daq.header.notes = 'Created using a GPOPS-II';


%set up Aux data
% auxdata.track = track ;
auxdata.vehicle = vehicle;
% auxdata.cost = cost;
%-------------------------------------------------------------------------%
%----------------------- Setup for Problem Bounds ------------------------%
%-------------------------------------------------------------------------%
bounds.phase.initialtime.lower  = s0; 
bounds.phase.initialtime.upper  = s0;
bounds.phase.finaltime.lower    = sf; 
bounds.phase.finaltime.upper    = sf;
bounds.phase.initialstate.lower = x0;
bounds.phase.initialstate.upper = x0;
bounds.phase.state.lower        = bounds.lbX; 
bounds.phase.state.upper        = bounds.ubX; 
bounds.phase.finalstate.lower   = bounds.lbX; 
bounds.phase.finalstate.upper   = bounds.ubX; 
bounds.phase.control.lower      = bounds.lbU; 
bounds.phase.control.upper      = bounds.ubU; 
bounds.phase.path.lower         = bounds.pathLower;
bounds.phase.path.upper         = bounds.pathUpper;


%-------------------------------------------------------------------------%
%---------------------- Provide Guess of Solution ------------------------%
%-------------------------------------------------------------------------%
guess.phase.time    = timeGuess; 
guess.phase.state   = stateGuess;
guess.phase.control = controlGuess;
% guess.phase.integral = trapz(timeGuess,1./stateGuess(:,1));



%-------------------------------------------------------------------------%
%----------Provide Mesh Refinement Method and Initial Mesh ---------------%
%-------------------------------------------------------------------------%
mesh = gpopsOptions.mesh;

%-------------------------------------------------------------------------%
%------------- Assemble Information into Problem Structure ---------------%        
%-------------------------------------------------------------------------%
setup.name                        = gpopsOptions.setup.name;
setup.nlp.solver                  = gpopsOptions.setup.nlp.solver;
setup.derivatives.supplier        = gpopsOptions.setup.derivatives.supplier;%'adigator';%'sparseFD'; %'adigator';%
setup.derivatives.derivativelevel = gpopsOptions.setup.derivatives.derivativelevel;
setup.scales.method               = gpopsOptions.setup.scales.method;
setup.method                      = gpopsOptions.setup.method;
setup.displaylevel                = gpopsOptions.setup.displaylevel;

setup.functions.continuous        = @fourwheelContinuous;
setup.functions.endpoint          = @fourwheelEndpoint;

setup.auxdata                     = auxdata;
setup.bounds                      = bounds;
setup.guess                       = guess;
setup.mesh                        = mesh; 

%% Adigator links:
% if strcmp(linkAdigatorFiles,'manual')
    adigatorfilenames = adigatorGenFiles4gpops2(setup);
%     return
% % end
% 
% switch linkAdigatorFiles
%     case 'manual'
        setup.adigatorgrd.continuous    = @fourwheelContinuousADiGatorGrd;
        setup.adigatorgrd.endpoint      = @fourwheelEndpointADiGatorGrd;
        setup.adigatorhes.continuous    = @fourwheelContinuousADiGatorHes;
        setup.adigatorhes.endpoint      = @fourwheelEndpointADiGatorHes;
% end




%% SOVLE
tic
output   = gpops2(setup);
toc

saveSnapshotofShortSeg = 'snapshot';
%%Deal w/ saving file
if ~isempty(saveSnapshotofShortSeg)
    save(saveSnapshotofShortSeg)
    ipoptfn = sprintf('%s_IPOPT.txt',strrep(saveSnapshotofShortSeg,'.mat',''));
    sysCommand = sprintf('mv quadCarIPOPTinfo.txt %s',ipoptfn);
    system(sysCommand);
end

%Post process
if ismember(output.result.nlpinfo,gpopsOptions.acceptableNlpOutpus)
    convergence = true;
    if ~isa(gpopsOptions.specifyMeshIterationSolution,'char')
        output.result = output.meshhistory(gpopsOptions.specifyMeshIterationSolution).result;
    end
else
    daq = [];
    convergence = false; 
    return
end

%Daq the solution
                            %    vx   vy   r  t
daqGpops = parseGpops2toDaq(output,stateNames,...
                              controlNames,...
                              indepVarName,...
                              units,...
                              names);
tempHeader = daq.header;                          
daq = catstruct(daqGpops,daq);
daq.gpopsSetup = setup;
daq.header = tempHeader;

%Calc algebraic states
daq = calculateAlgebraicStates(daq);


%%Deal w/ saving file - need to resave after editing
if ~isempty(saveSnapshotofShortSeg)
    save(saveSnapshotofShortSeg)
end

        

