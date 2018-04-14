c
global DAQ
loadFiles('driverB')
load('2015_Corvette_C7R')
DAQ{1}.vehicle = vehicle;
    
%%
fixDaqRemoveNanDataFromDaqFile
flipSignOnDaqChannels('ay');
%Deal with signs and other issues
chToFlipSign = {    'omegaWheel_L1'
    'omegaWheel_L2'
    'omegaWheel_R1'
    'omegaWheel_R2'
    'yawRate'};
flipSignOnDaqChannels(chToFlipSign)
daqUnitConversion('omegaWheel_L1',myConstants.radps2rpm,'rad/s')
daqUnitConversion('omegaWheel_R1',myConstants.radps2rpm,'rad/s')
daqUnitConversion('omegaWheel_L2',myConstants.radps2rpm,'rad/s')
daqUnitConversion('omegaWheel_R2',myConstants.radps2rpm,'rad/s')

%Design a filter
getChannelDataFromDaqFile(DAQ{1},{'time';'distance'});
dt = time(2)-time(1);
fs = 1/dt;
fn = fs/2;
ffilt = 5/fn;
[b,a] = cheby2(10,80,ffilt);

%Add steering
getChannelDataFromDaqFile(DAQ{1},'wheelSteerAngle','wheelPositionSuffix','standard');
delta = -1*mean([wheelSteerAngle_L1 wheelSteerAngle_R1],2)*myConstants.deg2rad;
DAQ{1}.rawData = addMathChannelsThatAreStandardChannels(DAQ{1}.rawData,'delta',delta,'','rad');
getChannelDataFromDaqFile(DAQ{1},'delta','applyFilter',{b,a})
DAQ{1}.rawData.delta.meas = delta;

[xData, yData] = prepareCurveData( time, delta );
ft = fittype( 'smoothingspline' );
opts = fitoptions( ft );
opts.SmoothingParam = 0.9;
[hFit, gof] = fit( xData, yData, ft, opts );   

plotDaqChannels('time','delta')
hold all
plot(time,ppval(hFit.p,time),'Color',plotColor(2))
swaSpline = hFit.p;
DAQ{1}.rawData.delta.meas = ppval(hFit.p,time);
DAQ{1}.rawData.wheelSteerAngle_L1 = DAQ{1}.rawData.delta;
DAQ{1}.rawData.wheelSteerAngle_R1 = DAQ{1}.rawData.delta;
% dDelta_dt = fnder(swaSpline,1);
% u1 = ppval(dDelta_dt,xData);
% DAQ{1}.rawData = addMathChannelsThatAreStandardChannels(DAQ{1}.rawData,'u1',u1,sprintf('Data filtered with %fHz cheby2. then smoothing spline fit at p=%f, then differentiated to form u1',ffilt*fn,opts.SmoothingParam),'rad/2');
% plotDaqChannels('distance','u1')

%%

mathChannelCalculateSSForcesAndMomentrsFromVehicleData([],...
    'vxConversion',1,...
    'axConversion',1,...
    'ayConversion',1,...
    'yawRateConversion',1,...
    'wheelAngleConversion',myConstants.deg2rad);


plotDaqChannelsAtEachWheelPosition('distance','fz')
%% Slip quanities
    mathChannelReffFromConstantValues([0.3424-0.026 0.3424-0.026 0.3557-0.01067 0.3557-0.01067])
mathChannelCalculateTireSlipAngleAndRatiosFromWheelAngle([],...
            'timeConversion',         1, ...                                    %Math needs [s]
            'speedConversion',        1, ...	                  %Math needs [m/s]
            'yawRateConversion',      1, ...	                   %Math needs [rad/s]
            'reffConversion',         1, ...                                 %Math needs [m]
            'deltaConversion' ,       1, ...
            'omegaConversion',        1);

% plotDaqChannelsAtEachWheelPosition('slipRatio','fx')        

%%
wheelPosSuffix = {'_L1','_R1','_L2','_R2'};
for iPos = 1:length(wheelPosSuffix)
    getChannelDataFromDaqFile([],{'kappa' 'slipRatio'},'wheelPositionSuffix',wheelPosSuffix{iPos},'applyFilter',{b,a},'returnVariablesWithWheelPositionSuffix',false);
    [xData, yData] = prepareCurveData( time, kappa );
    ft = fittype( 'smoothingspline' );
    opts = fitoptions( ft );
    opts.SmoothingParam = 0.6;
    [hFit, gof] = fit( xData, yData, ft, opts );   

    plotDaqChannels('time',['slipRatio' wheelPosSuffix{iPos}])
    hold all
    plot(time,ppval(hFit.p,time),'Color',plotColor(2))
    DAQ{1}.rawData.(['slipRatio' wheelPosSuffix{iPos}]).meas = ppval(hFit.p,time);

    
    
end

duplicateChannel('slipRatio_L1','kappa_L1')
duplicateChannel('slipRatio_R1','kappa_R1')
duplicateChannel('slipRatio_L2','kappa_L2')
duplicateChannel('slipRatio_R2','kappa_R2')
flipSignOnDaqChannels(...
    {'kappa_L1'
    'kappa_R1'
    'kappa_L2'
    'kappa_R2'});

%% Torque Rate estimate
% carFilename = '2015_Corvette_C7R.mat';
% load(carFilename);
% DAQ{1}.vehicle = vehicle;
% 
% getChannelDataFromDaqFile(DAQ{1},'fx','wheelPositionSuffix','standard','applyFilter',{b,a});
% plotDaqChannelsAtEachWheelPosition('distance','fx')
% reff_front = vehicle.tire_front.reff.meas;
% reff_rear = vehicle.tire_rear.reff.meas;
% torqueDemand = sum([fx_L1*reff_front fx_R1*reff_front fx_L2*reff_rear fx_R2*reff_rear],2);
% DAQ{1}.rawData = addMathChannelsThatAreStandardChannels(DAQ{1}.rawData,'torqueDemand',torqueDemand,'Calculated from estimated tire forces and model reff','N*m');
% 
% 
% [xData, yData] = prepareCurveData( time, torqueDemand );
% ft = fittype( 'smoothingspline' );
% opts = fitoptions( ft );
% opts.SmoothingParam = 0.95;
% [hFit, gof] = fit( xData, yData, ft, opts );   
% 
% plotDaqChannels('time','torqueDemand')
% hold all
% T = ppval(hFit.p,time);
% plot(time,T,'color',plotColor(2),'linestyle','--')
% legend('Data','fit')
% tSpline = hFit.p;
% dT_dt = fnder(tSpline,1);
% u2 = ppval(dT_dt,xData);
% DAQ{1}.rawData = addMathChannelsThatAreStandardChannels(DAQ{1}.rawData,'u2',u2,sprintf('Torque estimate Data filtered with %fHz cheby2. then smoothing spline fit at p=%f, then differentiated to form u1',ffilt*fn,opts.SmoothingParam),'N*m/s');
% 
% plotDaqChannels('distance','u2')
% 
[~,fn,~] = fileparts(DAQ{1}.header.filename);
DAQ{1}.header.filename = [fn '_EstInputs'];
saveFiles

