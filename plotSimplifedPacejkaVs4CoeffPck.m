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
fx_L1 = simplifiedPacejka(2000,0,kappa,vehicle.tire_front.coeff.meas);

hold all;plot(kappa,fx_L1);