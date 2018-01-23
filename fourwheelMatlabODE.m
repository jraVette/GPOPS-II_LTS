function dx = fourwheelMatlabODE(s,x,u,sSim,auxdata,varargin)

defaults = {'suppressOutput',false
            'tractionControlAbs',false
            'absSlipRatioFlag',0.13};
setDefaultsForVarargin(defaults,varargin);


u = interp1(sSim,u,s);

input.phase.time = s;
input.phase.state = x';
input.phase.control = u;
input.auxdata = auxdata;
phaseout = fourwheelContinuous(input);
dx = phaseout.dynamics';

%% Brute force ABS
kappa_L1 = phaseout.algebraicStates.slipRatio_L1.meas;
kappa_R1 = phaseout.algebraicStates.slipRatio_R1.meas;
kappa_L2 = phaseout.algebraicStates.slipRatio_L2.meas;
kappa_R2 = phaseout.algebraicStates.slipRatio_R2.meas;

absFlag = false;
if abs(kappa_L1)>absSlipRatioFlag || abs(kappa_R1)>absSlipRatioFlag || ...
   abs(kappa_L2)>absSlipRatioFlag || abs(kappa_R2)>absSlipRatioFlag
    input.phase.control = 0;
    phaseout = fourwheelContinuous(input);
    dx = phaseout.dynamics';
    absFlag = true;
end

if ~suppressOutput
    persistent count
    if isempty(count); count = 1; end
    
    %What vars do you want to print
    vars = {'s', 'x(1)'};
    headerText = [];
    dataText = [];
    for i = 1:length(vars)
        headerText = [headerText sprintf('%-5s ',vars{i})];
        dataText   = [dataText   sprintf('%5.1f ',eval(vars{i}))];
    end

    if count == 1
        fprintf('%s\n',headerText)
    end
    count = count+1;
    if count == 21
        count = 1;
    end
    
    if absFlag
        dataText = sprintf('%s ABS FLAG',dataText);
    end
    fprintf('%s\n',dataText)
    

end
