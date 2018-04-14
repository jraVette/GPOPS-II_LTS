function plotCostFunctionComponentsFromDaq(daq)
%Creation 23 Mar 2018 - Jeff Anderson

if ~exist('daq','var'); daq = []; end
DAQ = useLocalOrGlobalDaqFiles(daq);

for i = 1:length(DAQ)
    daq = DAQ{i};
    
    s   = writeDaqChannelsToMatrix(daq,'selectedChannels',daq.header.variableNames.indepVarName);
    x   = writeDaqChannelsToMatrix(daq,'selectedChannels',daq.header.variableNames.stateNames);
    u   = writeDaqChannelsToMatrix(daq,'selectedChannels',daq.header.variableNames.controlNames);
    
    w = daq.header.setup.auxdata.controlCost;
    uWeight = bsxfun(@times,u.^2,w);
    costTerms = trapz(s,uWeight);
    
    disp(' ');
    disp('_____________________________')
    disp('REGULARIZATION')
    disp(' ');
    for i = 1:length(costTerms)
        label{i} = sprintf('u_%i',i);
        fprintf('J_u%i, : %f\n',i,costTerms(i));
    end
    fprintf('Total regularization: %f\n',sum(costTerms));
    disp(' ');    
    disp('_____________________________')
    disp(' ');
    
    
    figure;
    subplot(1,2,1);
    
    pie(costTerms/sum(costTerms),label);
    title('Regularization Terms')
    
    subplot(1,2,2);

    tf = x(end,end);
    
    allCostComponents = [tf daq.header.setup.auxdata.regularizationCost*sum(costTerms)];
    pie(allCostComponents/sum(allCostComponents),{'Min Time','Regularization'});
    title('Cost Fcn')
    
    
        disp(' ');
    disp('_____________________________')
    disp('TOTAL COST')
    disp(' ');
    fprintf('Min Time Term: %f\n',tf);
    fprintf('Regularization: %f\n', allCostComponents(end))
    fprintf('Total Cost: %f\n', sum(allCostComponents))
    disp(' ');    
    disp('_____________________________')
    disp(' ');    
    


    
end