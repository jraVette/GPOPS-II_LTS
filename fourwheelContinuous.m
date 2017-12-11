function phaseout = fourwheelContinuous(input)


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
reff_front          = vehicle.tire_front.reff.meas;                        %Effective rolling radius front axle            [m]  
reff_rear           = vehicle.tire_rear.reff.meas;                        %Effective rolling radius front axle            [m]  
torqueBrakingRear   = vehicle.parameter.torqueDistBrakingRear.meas;        %Torque distrubution going to rear under braking [-] \in [0,1]
torqueDrivingRear   = vehicle.parameter.torqueDistDrivingRear.meas;        %Torque distrubution going to rear under driving [-] \in [0,1]
kd                  = vehicle.parameter.differentialFrictionCoeff.meas;    %Differential friction coeff                    [N*m*s/rad]


coeffFront          = vehicle.tire_front.coeff.meas;
coeffRear           = vehicle.tire_rear.coeff.meas; %STOPPED HERE

%% Constants that should be overridden
fz = 2000;
delta = 0;
Fax = 0;
k = 0;


%IndepVar
s                   = input.phase.time;

%States
vx                  = input.phase.state(:,1);
vy                  = input.phase.state(:,2);
r                   = input.phase.state(:,3);
omega_L1            = input.phase.state(:,4);
omega_R1            = input.phase.state(:,5);
omega_L2            = input.phase.state(:,6);
omega_R2            = input.phase.state(:,7);
T                   = input.phase.state(:,8);
ey                  = input.phase.state(:,9);
ePsi                = input.phase.state(:,10);

%Control
u2                  = input.phase.control(:,1)*5000;


%% Torque allocation - Power Train (based on the work in Tremlett)
tPlus   = 0.5+0.5*sin(atan(100*T));
tMinus  = 0.5-0.5*sin(atan(100*T));
kt = tPlus*torqueDrivingRear + tMinus*torqueBrakingRear;

T_drive_L1 = (1-kt).*(T)/(2);
T_drive_R1 = (1-kt).*(T)/(2);
T_drive_L2 = (kt).*(T)/(2) + kd*(omega_L2 - omega_R2);
T_drive_R2 = (kt).*(T)/(2) - kd*(omega_L2 - omega_R2);




%% Slip angle and ratios
kappa_L1    = -(1 + reff_front*-omega_L1./vx); %Note negative on slip ratio required because slip ratio defined like SAE to match Limebeer
kappa_R1    = -(1 + reff_front*-omega_R1./vx);
kappa_L2    = -(1 + reff_rear* -omega_L2./vx); 
kappa_R2    = -(1 + reff_rear* -omega_R2./vx);

[fx_L1, fy_L1] = simplifiedPacejka(fz,0,kappa_L1,coeffFront); %IFIX still need coeffRear to be actual rear tire
[fx_R1, fy_R1] = simplifiedPacejka(fz,0,kappa_R1,coeffFront);
[fx_L2, fy_L2] = simplifiedPacejka(fz,0,kappa_L2,coeffRear);
[fx_R2, fy_R2] = simplifiedPacejka(fz,0,kappa_R2,coeffRear);


%Sum forces and moment
FX = cos(delta).*(fx_L1 + fx_R1) - sin(delta).*(fy_L1 + fy_R1) + fx_L2 + fx_R2 + Fax;
FY = cos(delta).*(fy_L1 + fy_R1) + sin(delta).*(fx_L1 + fx_R1) + fy_L2 + fy_R2;


dvx_dt =  vy.*r + FX/m;
dvy_dt = -vx.*r + FY/m;
dr_dt = (a*(cos(delta).*(fy_R1 + fy_L1) + sin(delta).*(fx_R1 + fx_L1)) + ...
         wf*(fy_R1.*sin(delta) - fx_R1.*cos(delta)) - wr*fx_R2 + ...
         wf*(fx_L1.*cos(delta) - fy_L1.*sin(delta))  + wr*fx_L2 - b*(fy_R2 + fy_L2))/Izz;


%Wheel dynamics
domega_L2_dt = (T_drive_L2 - reff_front*fx_L2)./Jr_f; %IFIX still need J_r
domega_R2_dt = (T_drive_R2 - reff_front*fx_R2)./Jr_f;
domega_L1_dt = tMinus.*(T_drive_L1 - reff_rear*fx_L1)/Jr_f + tPlus.*(dvx_dt/reff_front); %Algebraic when free wheeling
domega_R1_dt = tMinus.*(T_drive_R1 - reff_rear*fx_R1)/Jr_f + tPlus.*(dvx_dt/reff_front);


% road dynamics dynamics
sDot = (vx.*cos(ePsi) - vy.*sin(ePsi))./(1-ey.*k);
dePsi_dt = r - k.*sDot;
dey_dt   = vx.*sin(ePsi) + vy.*cos(ePsi);



%% Outputs:
phaseout.dynamics = [dvx_dt,...
                     dvy_dt,...
                     dr_dt,...
                     domega_L1_dt,...
                     domega_R1_dt,...
                     domega_L2_dt,...
                     domega_R2_dt u2,...
                     dey_dt,...
                     dePsi_dt];
                 
                 
phaseout.integrand = 0.001*u2.^2;


phaseout.path = [kappa_L1,...
                 kappa_R1,...
                 kappa_L2,...
                 kappa_R2];
                

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

phaseout.algebraicStates.FX.meas = FX;
phaseout.algebraicStates.FY.meas = FY;

phaseout.algebraicStates.T_drive_L1.meas = T_drive_L1;
phaseout.algebraicStates.T_drive_R1.meas = T_drive_R1;
phaseout.algebraicStates.T_drive_L2.meas = T_drive_L2;
phaseout.algebraicStates.T_drive_R2.meas = T_drive_R2;


