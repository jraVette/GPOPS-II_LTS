kappa = -1:0.001:1; 

B = 10;
C = 2.75;
D = 6000;
E = 0.965;
fx_L1 = D*sin(C*  atan( B*kappa - E*(B*kappa - atan(B*kappa)))  );

figure;plot(kappa,fx_L1);


vehicleDirectory = fullfile(jatecPath,'Resources/In house code/Vehicle Parameters/');
carFilename = '2015_Corvette_C7R.mat';
fullVehicleFile = fullfile(vehicleDirectory,'Corvette',carFilename);
% carFileName = 'LimebeerF1Car.mat';
% fullVehicleFile = fullfile(vehicleDirectory,'Optimal Control Research',carFileName);
load(fullVehicleFile);
[fx_L1, Fy, muX, muY, FxMax, FyMax, kappa_n] = simplifiedPacejka(2000,0,kappa,vehicle.tire_front.coeff.meas);

hold all;plot(kappa,fx_L1);


%% Simplifed pck to 4 parameter
Fz = 2000;
muX_max = 2.96768625411826;
kappaMax = 0.11;
Qx = 1.9;

kappa_n = kappa./kappaMax;
rho = kappa_n;
Sx = pi/(2*atan(Qx));


muX = muX_max.*sin(coeff.Qx*atan(Sx*rho));

Fx = muX.*Fz.*kappa_n./(rho+eps)

hold all;
plot(kappa,Fx)


%Put it in 4 param method
FX = muX_max.*Fz.*kappa_n/(rho+eps).*sin(coeff.Qx*atan(Sx*rho)) 




hold all;
plot(kappa,FX)


figure;plot(kappa_n,FX)
