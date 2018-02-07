function phaseout = fourwheelContinuous(input)
%% Stats and controls, scale them appropiately
scaling = input.auxdata.scaling;
%IndepVar
s                   = input.phase.time./scaling.length;     

%States
vx                  = input.phase.state(:,1)./scaling.state(1);
vy                  = input.phase.state(:,2)./scaling.state(2);
r                   = input.phase.state(:,3)./scaling.state(3);
omega_L1            = input.phase.state(:,4)./scaling.state(4);
omega_R1            = input.phase.state(:,5)./scaling.state(5);
omega_L2            = input.phase.state(:,6)./scaling.state(6);
omega_R2            = input.phase.state(:,7)./scaling.state(7);
T                   = input.phase.state(:,8)./scaling.state(8);
ey                  = input.phase.state(:,9)./scaling.state(9);
ePsi                = input.phase.state(:,10)./scaling.state(10);
delta               = input.phase.state(:,11)./scaling.state(11);
%IF YOU ADD STATES, expand integral 


%Control
u1                  = input.phase.control(:,1)./scaling.control(1);
u2                  = input.phase.control(:,2)./scaling.control(2);


%Parameters
vehicle             = input.auxdata.vehicle;
m                   = vehicle.parameter.mass.meas;                         %Total vehicle mass                             [kg] 
Izz                 = vehicle.parameter.yawInertia.meas;                   %Yaw Inertia                                    [kg*m^2]
a                   = vehicle.parameter.a.meas ;                           %Length of front axle to CG                     [m]
b                   = vehicle.parameter.b.meas;                            %Length of rear axle to CG                      [m]
wf                  = vehicle.parameter.trackWidth_front.meas/2;           %1/2 of front track                             [m]
wr                  = vehicle.parameter.trackWidth_rear.meas/2;            %1/2 of rear track                              [m]
Jtire               = vehicle.tire_front.tireInertia.meas;     
Jwheel              = vehicle.tire_front.wheelInertia.meas;
Jr_f                = Jtire + Jwheel;                                      %Wheel inertia front wheel+tire                 [kg*m^2]
Jwheel              = vehicle.tire_rear.wheelInertia.meas;
Jr_r                = Jtire + Jwheel;  
reff_front          = vehicle.tire_front.reff.meas;                        %Effective rolling radius front axle            [m]  
reff_rear           = vehicle.tire_rear.reff.meas;                        %Effective rolling radius front axle            [m]  
torqueBrakingRear   = vehicle.parameter.torqueDistBrakingRear.meas;        %Torque distrubution going to rear under braking [-] \in [0,1]
torqueDrivingRear   = vehicle.parameter.torqueDistDrivingRear.meas;        %Torque distrubution going to rear under driving [-] \in [0,1]
kd                  = vehicle.parameter.differentialFrictionCoeff.meas;    %Differential friction coeff                    [N*m*s/rad]
CD                  = vehicle.parameter.coeffDrag.meas;                    %Coefficient of drag
CL                  = vehicle.parameter.coeffLift.meas;                    %Coefficient of lift
frontalArea         = vehicle.parameter.frontalArea.meas;                  %Frontal area                                   [m^2]
rho                 = vehicle.parameter.airDensity.meas;                   %Air Density                                    [kg/m^3]
a_a                 = vehicle.parameter.a_a.meas;                          %Dist to CP to front axle                       [m]
maxEnginePower      = vehicle.parameter.enginePower.meas;%*745.7010; %Max engine power                                [W]


coeffFront          = vehicle.tire_front.coeff.meas;
coeffRear           = vehicle.tire_rear.coeff.meas; 

%% Constants that should be overridden
fz = 2000;


%% Track
track = input.auxdata.track;
k = interp1(track.distance.meas,track.curvature.meas,s,'spline','extrap');

%% Aero
%Aero loads      
Faz =  0.5*CL*rho*frontalArea*vx.^2;
Fax = -0.5*CD*rho*frontalArea*vx.^2;


%% Torque allocation - Power Train (based on the work in Tremlett)
% tSign = T./abs(T);
torqueAllocationSlope  = input.auxdata.torqueAllocationSlope;
% tPlus   = 0.5+0.5*sin(atan(torqueAllocationSlope*tSign));
% tMinus  = 0.5-0.5*sin(atan(torqueAllocationSlope*tSign));
% tPlus = 0.5+05*tanh(torqueAllocationSlope*T);
% tMinus= 0.5-05*tanh(torqueAllocationSlope*T);
tPlus  = 0.5+0.5*(T./(sqrt(T.^2 + (1/torqueAllocationSlope)^2)));
% tMinus = 0.5-0.5*(T./(sqrt(T.^2 + (1/torqueAllocationSlope)^2)));

kt = tPlus*torqueDrivingRear + (1-tPlus)*torqueBrakingRear;

T_drive_L1 = (1-kt).*(T)/(2);
T_drive_R1 = (1-kt).*(T)/(2);
T_drive_L2 = (kt).*(T)/(2) + kd*(omega_L2 - omega_R2);
T_drive_R2 = (kt).*(T)/(2) - kd*(omega_L2 - omega_R2);




%% Slip angle and ratios
kappa_L1    = -(1 + reff_front*-omega_L1./vx); %Note negative on slip ratio required because slip ratio defined like SAE to match Limebeer
kappa_R1    = -(1 + reff_front*-omega_R1./vx);
kappa_L2    = -(1 + reff_rear* -omega_L2./vx); 
kappa_R2    = -(1 + reff_rear* -omega_R2./vx);

%Slip angles:
alpha_L1    = atan((-sin(delta).*(vx + r*wf) + cos(delta).*(vy + r*a))./(cos(delta).*(vx + r*wf) + sin(delta).*(vy + r*a)));
alpha_R1    = atan(( sin(delta).*(r*wf - vx) + cos(delta).*(vy + r*a))./(cos(delta).*(vx - r*wf) + sin(delta).*(vy + r*a)));
alpha_L2    = atan((vy - r*b)./(vx + r*wr));
alpha_R2    = atan((vy - r*b)./(vx - r*wr));


% [fx_L1, fy_L1] = simplifiedPacejka(fz,-alpha_L1,kappa_L1,coeffFront); 
% [fx_R1, fy_R1] = simplifiedPacejka(fz,-alpha_R1,kappa_R1,coeffFront);
% [fx_L2, fy_L2] = simplifiedPacejka(fz,-alpha_L2,kappa_L2,coeffRear);
% [fx_R2, fy_R2] = simplifiedPacejka(fz,-alpha_R2,kappa_R2,coeffRear);
% [fx_L1, fy_L1, muX_L1, muY_L1, FxMax_L1, FyMax_L1, kappa_n_L1, alpha_n_L1, rho_L1, eff_L1] = simplifiedPacejka(fz,-0,kappa_L1,coeffFront); 
% [fx_R1, fy_R1, muX_R1, muY_R1, FxMax_R1, FyMax_R1, kappa_n_R1, alpha_n_R1, rho_R1, eff_R1] = simplifiedPacejka(fz,-0,kappa_R1,coeffFront);
% [fx_L2, fy_L2, muX_L2, muY_L2, FxMax_L2, FyMax_L2, kappa_n_L2, alpha_n_L2, rho_L2, eff_L2] = simplifiedPacejka(fz,-0,kappa_L2,coeffRear);
% [fx_R2, fy_R2, muX_R2, muY_R2, FxMax_R2, FyMax_R2, kappa_n_R2, alpha_n_R2, rho_R2, eff_R2] = simplifiedPacejka(fz,-0,kappa_R2,coeffRear);
[fx_L1, fy_L1, muX_L1, muY_L1, ~, ~, kappa_n_L1, ~, rho_L1, eff_L1] = simplifiedPacejka(fz,-0,kappa_L1,coeffFront); 
[fx_R1, fy_R1, muX_R1, muY_R1, ~, ~, kappa_n_R1, ~, rho_R1, eff_R1] = simplifiedPacejka(fz,-0,kappa_R1,coeffFront);
[fx_L2, fy_L2, muX_L2, muY_L2, ~, ~, kappa_n_L2, ~, rho_L2, eff_L2] = simplifiedPacejka(fz,-0,kappa_L2,coeffRear);
[fx_R2, fy_R2, muX_R2, muY_R2, ~, ~, kappa_n_R2, ~, rho_R2, eff_R2] = simplifiedPacejka(fz,-0,kappa_R2,coeffRear);




%Sum forces and moment
FX = cos(delta).*(fx_L1 + fx_R1) - sin(delta).*(fy_L1 + fy_R1) + fx_L2 + fx_R2 + Fax;
FY = cos(delta).*(fy_L1 + fy_R1) + sin(delta).*(fx_L1 + fx_R1) + fy_L2 + fy_R2;


dvx_dt =  vy.*r + FX/m;
dvy_dt = -vx.*r + FY/m;
dr_dt = (a*(cos(delta).*(fy_R1 + fy_L1) + sin(delta).*(fx_R1 + fx_L1)) + ...
         wf*(fy_R1.*sin(delta) - fx_R1.*cos(delta)) - wr*fx_R2 + ...
         wf*(fx_L1.*cos(delta) - fy_L1.*sin(delta))  + wr*fx_L2 - b*(fy_R2 + fy_L2))/Izz;


%Wheel dynamics
domega_L2_dt = (T_drive_L2 - reff_front*fx_L2)./Jr_r; 
domega_R2_dt = (T_drive_R2 - reff_front*fx_R2)./Jr_r;

domega_L1_dt = (1-tPlus).*(T_drive_L1 - reff_rear*fx_L1)/Jr_f + tPlus.*(dvx_dt/reff_front); %Algebraic when free wheeling
domega_R1_dt = (1-tPlus).*(T_drive_R1 - reff_rear*fx_R1)/Jr_f + tPlus.*(dvx_dt/reff_front);


% road dynamics dynamics
sDot = (vx.*cos(ePsi) - vy.*sin(ePsi))./(1-ey.*k);
dePsi_dt = r - k.*sDot;
dey_dt   = vx.*sin(ePsi) + vy.*cos(ePsi);

%Engine
wheelSpeed = (omega_L2+omega_R2)/2; %Need positive number
percentEnginePowerUsed = ((T).*wheelSpeed)/maxEnginePower;
                 %Torque [N*m]*[rad/s]       [W]    
                 %Power [N*m/s = W]


%% Outputs:
phaseout.dynamics = [(dvx_dt./sDot).*(scaling.velocity/scaling.time/scaling.velocity),...
                     (dvy_dt./sDot).*(scaling.velocity/scaling.time/scaling.velocity),...
                     (dr_dt./sDot).*(scaling.angularVelocity/scaling.time/scaling.velocity),...
                     (domega_L1_dt./sDot).*(scaling.angularVelocity/scaling.time/scaling.velocity),...
                     (domega_R1_dt./sDot).*(scaling.angularVelocity/scaling.time/scaling.velocity),...
                     (domega_L2_dt./sDot).*(scaling.angularVelocity/scaling.time/scaling.velocity),...
                     (domega_R2_dt./sDot).*(scaling.angularVelocity/scaling.time/scaling.velocity),...
                     (u2./sDot).*(scaling.torque/scaling.time/scaling.velocity),...
                     (dey_dt./sDot).*(scaling.length/scaling.time/scaling.velocity),...
                     (dePsi_dt./sDot).*(scaling.angle/scaling.time/scaling.velocity),...
                     (u1./sDot).*(scaling.angle/scaling.time/scaling.velocity)];              
                 
% phaseout.integrand = 100*(1./sDot) + input.auxdata.controlWeight*u2.^2 + ... %worked with 100*minTime
%                      (vx-100).^2;%;*1e-1;
phaseout.integrand = (vx-100).^2 + vy.^2 + r.^2 + ey.^2 + ePsi.^2 + 1./sDot + input.auxdata.controlWeight*u2.^2;
timePenality = 1./sDot;
% stageCost = bsxfun(@times,bsxfun(@minus,input.phase.state,input.auxdata.stageCost.targetState).^2,input.auxdata.stageCost.weight);
controlCost = input.auxdata.controlWeight*u2.^2;

stageCost= (input.phase.state(:,1) - input.auxdata.stageCost.targetState(1)).^2*input.auxdata.stageCost.scaling(1) + ...
           (input.phase.state(:,2) - input.auxdata.stageCost.targetState(2)).^2*input.auxdata.stageCost.scaling(2) + ...
           (input.phase.state(:,3) - input.auxdata.stageCost.targetState(3)).^2*input.auxdata.stageCost.scaling(3) + ...
           (input.phase.state(:,4) - input.auxdata.stageCost.targetState(4)).^2*input.auxdata.stageCost.scaling(4) + ...
           (input.phase.state(:,5) - input.auxdata.stageCost.targetState(5)).^2*input.auxdata.stageCost.scaling(5) + ...
           (input.phase.state(:,6) - input.auxdata.stageCost.targetState(6)).^2*input.auxdata.stageCost.scaling(6) + ...
           (input.phase.state(:,7) - input.auxdata.stageCost.targetState(7)).^2*input.auxdata.stageCost.scaling(7) + ...
           (input.phase.state(:,8) - input.auxdata.stageCost.targetState(8)).^2*input.auxdata.stageCost.scaling(8) + ...
           (input.phase.state(:,9) - input.auxdata.stageCost.targetState(9)).^2*input.auxdata.stageCost.scaling(9) + ...
           (input.phase.state(:,10) - input.auxdata.stageCost.targetState(10)).^2*input.auxdata.stageCost.scaling(10) + ...
           (input.phase.state(:,11) - input.auxdata.stageCost.targetState(11)).^2*input.auxdata.stageCost.scaling(11);

phaseout.integral = timePenality*input.auxdata.minTimeCost + ...
                    controlCost + ...
                    stageCost*input.auxdata.stageCost.weight;
% phaseout.path = [kappa_L1,...
%                  kappa_R1,...
%                  kappa_L2,...
%                  kappa_R2];%;,...
%                  percentEnginePowerUsed-1];

phaseout.path = [kappa_n_L1 kappa_n_R1 kappa_n_L2 kappa_n_R2 percentEnginePowerUsed];

phaseout.algebraicStates.slipRatio_L1.meas = kappa_L1;
phaseout.algebraicStates.slipRatio_R1.meas = kappa_R1;
phaseout.algebraicStates.slipRatio_L2.meas = kappa_L2;
phaseout.algebraicStates.slipRatio_R2.meas = kappa_R2;

phaseout.algebraicStates.fx_L1.meas = fx_L1;
phaseout.algebraicStates.fx_R1.meas = fx_R1;
phaseout.algebraicStates.fx_L2.meas = fx_L2;
phaseout.algebraicStates.fx_R2.meas = fx_R2;

phaseout.algebraicStates.fy_L1.meas = fy_L1;
phaseout.algebraicStates.fy_R1.meas = fy_R1;
phaseout.algebraicStates.fy_L2.meas = fy_L2;
phaseout.algebraicStates.fy_R2.meas = fy_R2;

phaseout.algebraicStates.muX_L1.meas = muX_L1;
phaseout.algebraicStates.muX_R1.meas = muX_R1;
phaseout.algebraicStates.muX_L2.meas = muX_L2;
phaseout.algebraicStates.muX_R2.meas = muX_R2;

phaseout.algebraicStates.muY_L1.meas = muY_L1;
phaseout.algebraicStates.muY_R1.meas = muY_R1;
phaseout.algebraicStates.muY_L2.meas = muY_L2;
phaseout.algebraicStates.muY_R2.meas = muY_R2;

% phaseout.algebraicStates.FxMax_L1.meas = FxMax_L1;
% phaseout.algebraicStates.FxMax_R1.meas = FxMax_R1;
% phaseout.algebraicStates.FxMax_L2.meas = FxMax_L2;
% phaseout.algebraicStates.FxMax_R2.meas = FxMax_R2;
% 
% phaseout.algebraicStates.FyMax_L1.meas = FyMax_L1;
% phaseout.algebraicStates.FyMax_R1.meas = FyMax_R1;
% phaseout.algebraicStates.FyMax_L2.meas = FyMax_L2;
% phaseout.algebraicStates.FyMax_R2.meas = FyMax_R2;

phaseout.algebraicStates.kappa_n_L1.meas = kappa_n_L1;
phaseout.algebraicStates.kappa_n_R1.meas = kappa_n_R1;
phaseout.algebraicStates.kappa_n_L2.meas = kappa_n_L2;
phaseout.algebraicStates.kappa_n_R2.meas = kappa_n_R2;

% phaseout.algebraicStates.alpha_n_L1.meas = alpha_n_L1;
% phaseout.algebraicStates.alpha_n_R1.meas = alpha_n_R1;
% phaseout.algebraicStates.alpha_n_L2.meas = alpha_n_L2;
% phaseout.algebraicStates.alpha_n_R2.meas = alpha_n_R2;

phaseout.algebraicStates.rho_L1.meas = rho_L1;
phaseout.algebraicStates.rho_R1.meas = rho_R1;
phaseout.algebraicStates.rho_L2.meas = rho_L2;
phaseout.algebraicStates.rho_R2.meas = rho_R2;

phaseout.algebraicStates.eff_L1.meas = eff_L1;
phaseout.algebraicStates.eff_R1.meas = eff_R1;
phaseout.algebraicStates.eff_L2.meas = eff_L2;
phaseout.algebraicStates.eff_R2.meas = eff_R2;


phaseout.algebraicStates.FX.meas = FX;
phaseout.algebraicStates.FY.meas = FY;
phaseout.algebraicStates.ax.meas = FX./m;
phaseout.algebraicStates.ay.meas = FY./m;

phaseout.algebraicStates.T_drive_L1.meas = T_drive_L1;
phaseout.algebraicStates.T_drive_R1.meas = T_drive_R1;
phaseout.algebraicStates.T_drive_L2.meas = T_drive_L2;
phaseout.algebraicStates.T_drive_R2.meas = T_drive_R2;

phaseout.algebraicStates.kt.meas = kt;
phaseout.algebraicStates.diffTorqueTransfer.meas = kd*(omega_L2 - omega_R2);
phaseout.algebraicStates.tPlus.meas =tPlus;
phaseout.algebraicStates.tMinus.meas =(1-tPlus);

phaseout.algebraicStates.enginePercent.meas = percentEnginePowerUsed;

phaseout.algebraicStates.vx_unscaled.meas = vx;
phaseout.algebraicStates.vy_unscaled.meas = vy;
phaseout.algebraicStates.r_unscaled.meas = r;
phaseout.algebraicStates.omega_L1_unscaled.meas = omega_L1;
phaseout.algebraicStates.omega_R1_unscaled.meas = omega_R1;
phaseout.algebraicStates.omega_L2_unscaled.meas = omega_L2;
phaseout.algebraicStates.omega_R2_unscaled.meas = omega_R2;
phaseout.algebraicStates.T_unscaled.meas = T;
phaseout.algebraicStates.ey_unscaled.meas = ey;
phaseout.algebraicStates.ePsi_unscaled.meas = ePsi;
phaseout.algebraicStates.delta_unscaled.meas = delta;


