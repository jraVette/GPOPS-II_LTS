function daq = runDoubleTrackMatlabOde(s,x0,u,auxdata)
clc

if nargin == 0
    carFileName = 'LimebeerF1Car.mat';
    vehicleDirectory = fullfile(jatecPath,'Resources/In house code/Vehicle Parameters/');
    fullVehicleFile = fullfile(vehicleDirectory,'Optimal Control Research',carFileName);
    load(fullVehicleFile);
    s0          = 0;                                             %[m] s
    sf          = 10;
    s = s0:1:sf;
    u = 1000*ones(size(s));
    vx0 = 10;
    omega_front0 = -vx0./vehicle.tire_front.reff.meas;
    x0 = [vx0 omega_front0];
    auxdata.vehicle = vehicle;
end

%% Run simulation
u = rowVector(u);
sSpace = s;
sSpan = [s(1) s(end)]; %#ok s will be an array



options = odeset('MaxStep',1);
[sSol,xSol] = ode45(@(s,x)fourwheelMatlabODE(s,x,u,sSpace,auxdata), sSpan, x0,options);
% stateNames = { 'vx','omegaWheel_L1'};%,'torqueDemand'};
% controlNames = {'T'};
% units =      {'s';'m/s';'rad/s';'N*m'};
% names = {'Time';'Vx';'Wheel Speed Left Front';'Torque'};

stateNames = { 'vx'};%,'torqueDemand'};
controlNames = {'slipRatio'};
units =      {'s';'m/s';''};
names = {'Time';'Vx';'Slip Ratio'};
     
%Resample to requested input
SSOL = s;
for i = 1:size(xSol,2)
    XSOL(:,i) = interp1(sSol,xSol(:,i),SSOL);
end


daq = parseMatlabOdeOutput(SSOL,XSOL,u,'time',stateNames,controlNames,units,names);
gpopsOutput.result.solution.phase.time = SSOL;
gpopsOutput.result.solution.phase.state = XSOL;
gpopsOutput.result.solution.phase.control = u;
gpopsOutput = addNotesToDaqFile(gpopsOutput,'Just added this to mimic gpops output for calculating algebraic states');
daq.gpopsOutput = gpopsOutput; 
daq.vehicle = auxdata.vehicle;
if nargin == 0
    daq.vehicle = auxdata.vehicle;
%     daq = mathChannelCalculateAlgebraicStates(daq);

     quickPlotTable(daq)

end
