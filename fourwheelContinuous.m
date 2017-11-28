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
T_drive_L1          = input.phase.state(:,3);

%Control
u2                  = input.phase.control(:,1)*5000;
% T_drive_L1          = input.phase.control(:,1)*5000;
% kappa_L1          = input.phase.control(:,1);

%% Dynamic System
kappa_L1    = -(1 + reff_f*-omega_L1./vx);


% kappa_L1 = (-reff_f*omega_L1 - vx)./vx;
% fx_L1 = simplifiedPacejka(2000,0,kappa_L1,coeffFront);

% kappa_L1 = -1:0.0001:1


% B = 10;
% C = 2.75;
% D = 6000;
% E = 0.965;
% fx_L1 = D*sin(C*  atan( B*kappa_L1 - E*(B*kappa_L1 - atan(B*kappa_L1)))  );

%Simplified pck
Fz = 2000;
muX_max = 2.96768625411826;
kappaMax = 0.11;
Qx = 1.9;
kappa_n = kappa_L1./kappaMax;
rho = kappa_n;
Sx = pi/(2*atan(Qx));
muX = muX_max.*sin(Qx*atan(Sx*rho));
fx_L1 = muX.*Fz.*kappa_n./(rho+eps);


dvx_dt = fx_L1./m;

domega_L1_dt = (T_drive_L1 - reff_f*fx_L1)./Jr_f;


phaseout.dynamics = [dvx_dt domega_L1_dt u2];
% phaseout.dynamics = [dvx_dt ];
% phaseout.integrand = u2.^2;
phaseout.path = kappa_L1;

phaseout.algebraicStates.slipRatio_L1.meas = kappa_L1;
phaseout.algebraicStates.fx_L1.meas = fx_L1;


% disp(kappa_L1)
