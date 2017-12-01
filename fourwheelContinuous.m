function phaseout = fourwheelContinuous(input)


%Parameters
vehicle             = input.auxdata.vehicle;
m                   = vehicle.parameter.mass.meas;                         %Total vehicle mass                             [kg] 
Jtire               = vehicle.tire_front.tireInertia.meas;     
Jwheel              = vehicle.tire_front.wheelInertia.meas;
Jr_f                = Jtire + Jwheel;                                      %Wheel inertia front wheel+tire                 [kg*m^2]
reff_f              = vehicle.tire_front.reff.meas;                        %Effective rolling radius front axle            [m]  
reff_r              = vehicle.tire_rear.reff.meas;                        %Effective rolling radius front axle            [m]  

coeffFront          = vehicle.tire_front.coeff.meas;
coeffRear          = vehicle.tire_front.coeff.meas;

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

%% Torque allocation
T_drive_L1 = T/4;
T_drive_R1 = T/4;
T_drive_L2 = T/4;
T_drive_R2 = T/4;



%% Dynamic System
kappa_L1    = -(1 + reff_f*-omega_L1./vx);
kappa_R1    = -(1 + reff_f*-omega_R1./vx);
kappa_L2    = -(1 + reff_r*-omega_L2./vx);
kappa_R2    = -(1 + reff_r*-omega_R2./vx);

fx_L1 = simplifiedPacejka(2000,0,kappa_L1,coeffFront); %IFIX still need coeffRear to be actual rear tire
fx_R1 = simplifiedPacejka(2000,0,kappa_R1,coeffFront);
fx_L2 = simplifiedPacejka(2000,0,kappa_L2,coeffRear);
fx_R2 = simplifiedPacejka(2000,0,kappa_R2,coeffRear);

FX = fx_L1 + fx_R1 + fx_L2 + fx_R2;


dvx_dt = FX./m;

domega_L1_dt = (T_drive_L1 - reff_f*fx_L1)./Jr_f; %IFIX still need J_r
domega_R1_dt = (T_drive_R1 - reff_f*fx_R1)./Jr_f;
domega_L2_dt = (T_drive_L2 - reff_f*fx_L2)./Jr_f;
domega_R2_dt = (T_drive_R2 - reff_f*fx_R2)./Jr_f;


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


% disp(kappa_L1)
