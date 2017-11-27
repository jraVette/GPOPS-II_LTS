function phaseout = fourwheelContinuous(input)


%Parameters
vehicle             = input.auxdata.vehicle;
m                   = vehicle.parameter.mass.meas;                         %Total vehicle mass                             [kg] 
Jtire               = vehicle.tire_front.tireInertia.meas;     
Jwheel              = vehicle.tire_front.wheelInertia.meas;
Jr_f                = Jtire + Jwheel;                                      %Wheel inertia front wheel+tire                 [kg*m^2]
reff_f              = vehicle.tire_front.reff.meas;                        %Effective rolling radius front axle            [m]  

coeffFront          = vehicle.tire_front.coeff.meas;

%IndepVar
s                   = input.phase.time;

%States
vx                  = input.phase.state(:,1);
omega_L1            = input.phase.state(:,2);

%Control
T_drive_L1          = input.phase.control(:,1)*5000;
% kappa_L1          = input.phase.control(:,1);

%% Dynamic System
kappa_L1    = -(1 + reff_f*-omega_L1./vx);


% kappa_L1 = (-reff_f*omega_L1 - vx)./vx;
% fx_L1 = simplifiedPacejka(2000,0,kappa_L1,coeffFront);

% kappa_L1 = -1:0.0001:1
B = 10;
C = 2;
D = 2000;
E = 1;
fx_L1 = D*sin(C*  atan( B*kappa_L1 - E*(B*kappa_L1 - atan(B*kappa_L1)))  );


dvx_dt = fx_L1./m;

domega_L1_dt = (T_drive_L1 - reff_f*fx_L1)./Jr_f;


phaseout.dynamics = [dvx_dt domega_L1_dt];
% phaseout.dynamics = [dvx_dt ];
% phaseout.integrand = 1./vx;
phaseout.path = kappa_L1;

phaseout.algebraicStates.slipRatio_L1.meas = kappa_L1;
phaseout.algebraicStates.fx_L1.meas = fx_L1;


% disp(kappa_L1)
