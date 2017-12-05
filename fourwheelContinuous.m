function phaseout = fourwheelContinuous(input)


%Parameters
vehicle             = input.auxdata.vehicle;
m                   = vehicle.parameter.mass.meas;                         %Total vehicle mass                             [kg] 
Jtire               = vehicle.tire_front.tireInertia.meas;     
Jwheel              = vehicle.tire_front.wheelInertia.meas;
Jr_f                = Jtire + Jwheel;                                      %Wheel inertia front wheel+tire                 [kg*m^2]
reff_front          = vehicle.tire_front.reff.meas;                        %Effective rolling radius front axle            [m]  
reff_rear           = vehicle.tire_rear.reff.meas;                        %Effective rolling radius front axle            [m]  
torqueBrakingRear   = vehicle.parameter.torqueDistBrakingRear.meas;        %Torque distrubution going to rear under braking [-] \in [0,1]
torqueDrivingRear   = vehicle.parameter.torqueDistDrivingRear.meas;        %Torque distrubution going to rear under driving [-] \in [0,1]


coeffFront          = vehicle.tire_front.coeff.meas;
coeffRear           = vehicle.tire_rear.coeff.meas; %STOPPED HERE

%IndepVar
s                   = input.phase.time;

%States
vx                  = input.phase.state(:,1);
omega_L1            = input.phase.state(:,2);
omega_R1            = input.phase.state(:,3);
omega_L2            = input.phase.state(:,4);
omega_R2            = input.phase.state(:,5);
T                   = input.phase.state(:,6);

%Control
u2                  = input.phase.control(:,1)*5000;
% T_drive_L1          = input.phase.control(:,1)*5000;
% kappa_L1          = input.phase.control(:,1);

%% Torque allocation - Power Train (based on the work in Tremlett)
tPlus   = 0.5+0.5*sin(atan(100*T));
tMinus  = 0.5-0.5*sin(atan(100*T));
kt = tPlus*torqueDrivingRear + tMinus*torqueBrakingRear;

T_drive_L1 = (1-kt).*(T)/(2);
T_drive_R1 = (1-kt).*(T)/(2);
T_drive_L2 = (kt).*(T)/(2);% + kd*(omega_L2 - omega_R2);
T_drive_R2 = (kt).*(T)/(2);% - kd*(omega_L2 - omega_R2);



%% Dynamic System
kappa_L1    = -(1 + reff_front*-omega_L1./vx);
kappa_R1    = -(1 + reff_front*-omega_R1./vx);
kappa_L2    = -(1 + reff_rear*-omega_L2./vx); %fix reff
kappa_R2    = -(1 + reff_rear*-omega_R2./vx);

fx_L1 = simplifiedPacejka(2000,0,kappa_L1,coeffFront); %IFIX still need coeffRear to be actual rear tire
fx_R1 = simplifiedPacejka(2000,0,kappa_R1,coeffFront);
fx_L2 = simplifiedPacejka(2000,0,kappa_L2,coeffRear);
fx_R2 = simplifiedPacejka(2000,0,kappa_R2,coeffRear);

FX = fx_L1 + fx_R1 + fx_L2 + fx_R2;


dvx_dt = FX./m;

domega_L2_dt = (T_drive_L2 - reff_front*fx_L2)./Jr_f; %IFIX still need J_r
domega_R2_dt = (T_drive_R2 - reff_front*fx_R2)./Jr_f;
domega_L1_dt = tMinus.*(T_drive_L1 - reff_rear*fx_L1)/Jr_f + tPlus.*(dvx_dt/reff_front); %Algebraic when free wheeling
domega_R1_dt = tMinus.*(T_drive_R1 - reff_rear*fx_R1)/Jr_f + tPlus.*(dvx_dt/reff_front);


phaseout.dynamics = [dvx_dt domega_L1_dt domega_R1_dt domega_L2_dt domega_R2_dt u2];
% phaseout.dynamics = [dvx_dt ];
% phaseout.integrand = u2.^2;
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

phaseout.algebraicStates.T_drive_L1.meas = T_drive_L1;
phaseout.algebraicStates.T_drive_R1.meas = T_drive_R1;
phaseout.algebraicStates.T_drive_L2.meas = T_drive_L2;
phaseout.algebraicStates.T_drive_R2.meas = T_drive_R2;


% disp(kappa_L1)
