function phaseout = fourwheelContinuous(input)
%% Stats and controls, scale them appropiately
scaling = input.auxdata.scaling;
%IndepVar
s                   = input.phase.time./scaling.length;     

%States
vx                  = input.phase.state(:,1)./scaling.state(1);
vy                  = input.phase.state(:,2)./scaling.state(2);
r                   = input.phase.state(:,3)./scaling.state(3);
ey                  = input.phase.state(:,4)./scaling.state(4);
ePsi                = input.phase.state(:,5)./scaling.state(5);
t                   = input.phase.state(:,6)./scaling.state(6);


%Control
delta               = input.phase.control(:,1)./scaling.control(1);
kappa_L1            = input.phase.control(:,2)./scaling.control(2);
kappa_R1            = input.phase.control(:,3)./scaling.control(3);
kappa_L2            = input.phase.control(:,4)./scaling.control(4);
kappa_R2            = input.phase.control(:,5)./scaling.control(5);
fz_L1               = input.phase.control(:,6)./scaling.control(6);
fz_R1               = input.phase.control(:,7)./scaling.control(7);
fz_L2               = input.phase.control(:,8)./scaling.control(8);
fz_R2               = input.phase.control(:,9)./scaling.control(9);


%Parameters
g = 9.81;
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
D                   = vehicle.parameter.rollStiffnessDistribution.meas;    %Roll stiffness distribution of the front axle  
h                   = vehicle.parameter.hcg.meas;                          %Height of CG                                   [m]




maxEnginePower = input.auxdata.engMult*maxEnginePower;

coeffFront          = vehicle.tire_front.coeff.meas;
coeffFront.muX1 = input.auxdata.muMultX*coeffFront.muX1;
coeffFront.muX2 = input.auxdata.muMultX*coeffFront.muX2;
coeffFront.muY1 = input.auxdata.muMultY*coeffFront.muY1;
coeffFront.muY2 = input.auxdata.muMultY*coeffFront.muY2;


coeffRear           = vehicle.tire_rear.coeff.meas; 
coeffRear.muX1 = input.auxdata.muMultX*coeffRear.muX1;
coeffRear.muX2 = input.auxdata.muMultX*coeffRear.muX2;
coeffRear.muY1 = input.auxdata.muMultY*coeffRear.muY1;
coeffRear.muY2 = input.auxdata.muMultY*coeffRear.muY2;

%% Track
track = input.auxdata.track;
k = interp1(track.distance.meas,track.curvature.meas,s+input.auxdata.currentDistance,'spline','extrap');

%% Aero
%Aero loads      
Faz =  0.5*CL*rho*frontalArea*vx.^2;
Fax = -0.5*CD*rho*frontalArea*vx.^2;




%% Slip angle and ratios
%Slip Ratios 
% delta = -delta;
omega_L1 = -(-kappa_L1-1).*(cos(delta).*(vx + r*wf) + sin(delta).*(r*a + vy))/reff_front; 
omega_R1 = -(-kappa_R1-1).*(cos(delta).*(vx - r*wf) + sin(delta).*(r*a + vy))/reff_front;
omega_L2 = -(-kappa_L2-1).*(vx + r*wr)/reff_rear; 
omega_R2 = -(-kappa_R2-1).*(vx - r*wr)/reff_rear; 

%Slip angles:
alpha_L1    = atan((-sin(delta).*(vx + r*wf) + cos(delta).*(vy + r*a))./(cos(delta).*(vx + r*wf) + sin(delta).*(vy + r*a)));
alpha_R1    = atan(( sin(delta).*(r*wf - vx) + cos(delta).*(vy + r*a))./(cos(delta).*(vx - r*wf) + sin(delta).*(vy + r*a)));
alpha_L2    = atan((vy - r*b)./(vx + r*wr));
alpha_R2    = atan((vy - r*b)./(vx - r*wr));

[fx_L1, fy_L1, muX_L1, muY_L1, ~, ~, kappa_n_L1, ~, rho_L1, eff_L1] = simplifiedPacejka(fz_L1,alpha_L1,kappa_L1,coeffFront); 
[fx_R1, fy_R1, muX_R1, muY_R1, ~, ~, kappa_n_R1, ~, rho_R1, eff_R1] = simplifiedPacejka(fz_R1,alpha_R1,kappa_R1,coeffFront);
[fx_L2, fy_L2, muX_L2, muY_L2, ~, ~, kappa_n_L2, ~, rho_L2, eff_L2] = simplifiedPacejka(fz_L2,alpha_L2,kappa_L2,coeffRear);
[fx_R2, fy_R2, muX_R2, muY_R2, ~, ~, kappa_n_R2, ~, rho_R2, eff_R2] = simplifiedPacejka(fz_R2,alpha_R2,kappa_R2,coeffRear);

fy_L1 = fy_L1;
fy_R1 = fy_R1;
fy_L2 = fy_L2;
fy_R2 = fy_R2;

%Sum forces and moment
delta = -delta;
FX = cos(delta).*(fx_L1 + fx_R1) - sin(delta).*(fy_L1 + fy_R1) + fx_L2 + fx_R2 + Fax;
FY = cos(delta).*(fy_L1 + fy_R1) + sin(delta).*(fx_L1 + fx_R1) + fy_L2 + fy_R2;


dvx_dt =  vy.*r + FX/m;
dvy_dt = -vx.*r + FY/m;
dr_dt = (a*(cos(delta).*(fy_R1 + fy_L1) + sin(delta).*(fx_R1 + fx_L1)) + ...
         wf*(fy_R1.*sin(delta) - fx_R1.*cos(delta)) - wr*fx_R2 + ...
         wf*(fx_L1.*cos(delta) - fy_L1.*sin(delta))  + wr*fx_L2 - b*(fy_R2 + fy_L2))/Izz;



%% Road dynamics dynamics %NOTE I FLIPPED LATERAL SIGN
sDot = (vx.*cos(ePsi) - -vy.*sin(ePsi))./(1-ey.*k);
dePsi_dt = -r - k.*sDot;
dey_dt   = vx.*sin(ePsi) + -vy.*cos(ePsi);


%% Engine limits
T = fx_L2*reff_rear + fx_R2*reff_rear;
wheelSpeed = (omega_L2+omega_R2)/2; %Need positive number
percentEnginePowerUsed = ((T).*wheelSpeed)/maxEnginePower;


%% Outputs:
phaseout.dynamics = [(dvx_dt./sDot).*(scaling.velocity/scaling.time/scaling.velocity),...
                     (dvy_dt./sDot).*(scaling.velocity/scaling.time/scaling.velocity),...
                     (dr_dt./sDot).*(scaling.angularVelocity/scaling.time/scaling.velocity),...
                     (dey_dt./sDot).*(scaling.length/scaling.time/scaling.velocity),...
                     (dePsi_dt./sDot).*(scaling.angle/scaling.time/scaling.velocity),...
                     1./sDot/scaling.velocity];              
%  phaseout.integrand = 1000*(1./sDot);                
J1 =  input.auxdata.controlCost(1)*input.phase.control(:,1).^2 + ...
      input.auxdata.controlCost(2)*input.phase.control(:,2).^2 + ...
      input.auxdata.controlCost(3)*input.phase.control(:,3).^2 + ...
      input.auxdata.controlCost(4)*input.phase.control(:,4).^2 + ...
      input.auxdata.controlCost(5)*input.phase.control(:,5).^2 + ...
      input.auxdata.controlCost(6)*input.phase.control(:,6).^2 + ...
      input.auxdata.controlCost(7)*input.phase.control(:,7).^2 + ...
      input.auxdata.controlCost(8)*input.phase.control(:,8).^2 + ...
      input.auxdata.controlCost(9)*input.phase.control(:,9).^2;
J2 =  input.auxdata.controlCostVx(1)*input.phase.control(:,1).^2 + ...
      input.auxdata.controlCostVx(2)*input.phase.control(:,2).^2 + ...
      input.auxdata.controlCostVx(3)*input.phase.control(:,3).^2 + ...
      input.auxdata.controlCostVx(4)*input.phase.control(:,4).^2 + ...
      input.auxdata.controlCostVx(5)*input.phase.control(:,5).^2 + ...
      input.auxdata.controlCostVx(6)*input.phase.control(:,6).^2 + ...
      input.auxdata.controlCostVx(7)*input.phase.control(:,7).^2 + ...
      input.auxdata.controlCostVx(8)*input.phase.control(:,8).^2 + ...
      input.auxdata.controlCostVx(9)*input.phase.control(:,9).^2;
phaseout.integrand = [J1 J2];

% stateTracking = input.auxdata.Q(1)*(interp1(input.auxdata.refTime,input.auxdata.refStates(:,1),s+input.auxdata.currentDistance,'linear','extrap') - input.phase.state(:,1)./scaling.state(1)).^2 + ...
%                 input.auxdata.Q(2)*(interp1(input.auxdata.refTime,input.auxdata.refStates(:,2),s+input.auxdata.currentDistance,'linear','extrap') - input.phase.state(:,2)./scaling.state(2)).^2 + ...
%                 input.auxdata.Q(3)*(interp1(input.auxdata.refTime,input.auxdata.refStates(:,3),s+input.auxdata.currentDistance,'linear','extrap') - input.phase.state(:,3)./scaling.state(3)).^2 + ...
%                 input.auxdata.Q(4)*(interp1(input.auxdata.refTime,input.auxdata.refStates(:,4),s+input.auxdata.currentDistance,'linear','extrap') - input.phase.state(:,4)./scaling.state(4)).^2 + ...
%                 input.auxdata.Q(5)*(interp1(input.auxdata.refTime,input.auxdata.refStates(:,5),s+input.auxdata.currentDistance,'linear','extrap') - input.phase.state(:,5)./scaling.state(5)).^2;

% phaseout.integrand = controlPenalty+ stateTracking;

Tbrake_f = fx_L1*reff_front + fx_R1*reff_front;
Tbrake_r = fx_L2*reff_rear + fx_R2*reff_rear;

Tbrake_f = -(-Tbrake_f+sqrt(Tbrake_f.^2+1e-3))/2;
Tbrake_r = -(-Tbrake_r+sqrt(Tbrake_r.^2+1e-3))/2;


phaseout.path = [...
    fz_L1 + fz_R1 +  fz_L2  + fz_R2 + m*g + Faz, ...
     wr*(fz_L2 - fz_R2) + wf*(fz_L1 - fz_R1) + h*FY, ...
    b*(fz_R2 + fz_L2) - a*(fz_R1 + fz_L1) + h*FX + (a_a - a)*Faz, ...
    D*(fz_R1 + fz_R2 - fz_L1 - fz_L2) - (fz_R1 - fz_L1),...    
    (omega_R1+sqrt(omega_R1.^2+1e-3)/2).*(omega_L1+sqrt(omega_L1.^2+1e-3)/2).*(fx_R1 - fx_L1), ...
    -kd*(omega_L2 - omega_R2) - reff_rear*(fx_L2 - fx_R2),...
    percentEnginePowerUsed];
    
% (omega_R1+sqrt(omega_R1.^2+1e-3)/2).*(omega_L1+sqrt(omega_L1.^2+1e-3)/2).*(fx_R1 - fx_L1), ...
    %Tbrake_r./(Tbrake_f+Tbrake_r) - torqueBrakingRear, ...
%     percentEnginePowerUsed];


