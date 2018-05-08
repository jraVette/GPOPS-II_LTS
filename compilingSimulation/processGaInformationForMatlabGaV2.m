function [bestScore,bestPopulation,allScore,allPopulation,nonSortedScore,nonSortedPopulation,bestIterateInfo,allIterateInfo] = processGaInformationForMatlabGaV2(varargin)
%This funciton will use the outputfcn of the ga and save the structure
%gaInfo in a file.  The field gaInfo.generation(iState) will be the state
%variable directly from matlab (see ga output fcn for details on that
%structure).  We'll access the score and population information.
%INPTUS:
%    varargin - to set defaults
%OUTPUTS:
%    bestScore
%    bestPopulation
%    allScore
%    allPopulations
%
%Creation: 28 Sep 2016 - Jeff Anderson
%Updated:  17 Oct 2016 - Jeff Anderson - I started calling state generation
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


defaults = {'filename','gui'  %Used to specify the file or cell of files to get metrics for
            'suppressPlot',false
            'processGeneration',[]};%leave empty for all, or put integer list for the generations to include
setDefaultsForVarargin(defaults,varargin)

%GUI file selector
if strcmp(filename,'gui')
    [~,~,fullFilenames,justFilenames] = getFilenamesAndPath('gui','*.mat','dialogTitle','Choose files with gaSave variable');
%User specifed files    
else
    if isa(filename,'char')
        fullFilenames{1} = filename;
        justFilenames{1} = filename; %Used for plotting
    else
        fullFilenames = filename;
        justFilenames = filename; %Used for plotting
    end
end

%% Loop through files
for iFile = 1:length(fullFilenames)
    clear gaInfo
    load(fullFilenames{iFile})
    
    %Make sure this is a ga file
    if exist('gaInfo','var')
         nGen = length(gaInfo.generation);
    else
        error('gaInfo variable not detected');
    end
    
    
    %Initialize variables
    allScore = [];
    allPopulation = [];
    if isempty(processGeneration); processGeneration = 1:length(gaInfo.generation);end
        
    %Loop throught the GA's and grab the information
    for iGen = 1:nGen
        if ~isempty(gaInfo.generation(iGen).scores) && ismember(iGen,processGeneration)
            scoreGen{iGen} = gaInfo.generation(iGen).scores;
            allScore = [allScore; gaInfo.generation(iGen).scores];
            allPopulation = [allPopulation; gaInfo.generation(iGen).Population];
        end
    end
    
    nonSortedScore = allScore;
    nonSortedPopulation = allPopulation;
    
    [~,ind] = sort(allScore,'ascend');
    bestPopulation = allPopulation(ind(1),:);
    bestScore = allScore(ind(1));
    allPopulation = allPopulation(ind,:);
    allScore = allScore(ind);
    
    %Find the best iterate in the gaInfo file
    count = 1; %Conter for the best iterates
    for iGen = 1:nGen
        if ~isempty(gaInfo.generation(iGen).scores)
            genScore = gaInfo.generation(iGen).scores;
            genPop = gaInfo.generation(iGen).Population;
            if ismember(bestScore,genScore)
                [~,ind] = ismember(bestScore,genScore,'rows');
                bestIterateInfo(count).generation = iGen;
                bestIterateInfo(count).iterate = ind;
                bestIterateInfo(count).score = genScore(ind);
                bestIterateInfo(count).population = genPop(ind,:);
                count = count+1;
            end
        end
    end    
    
    %Get info for all other scores
    uniqueScores = unique(allScore);
    for iScore = 1:length(uniqueScores)
        iterateInfo.generation = [];
        iterateInfo.iterate    = [];
        iterateInfo.score      = [];
        iterateInfo.population = [];
        for iGen = 1:nGen
            if ~isempty(gaInfo.generation(iGen).scores)
                genScore = gaInfo.generation(iGen).scores;
                genPop = gaInfo.generation(iGen).Population;
                if ismember(uniqueScores(iScore),genScore)
                    [~,ind] = ismember(uniqueScores(iScore),genScore,'rows');
                    iterateInfo.generation    = [iterateInfo.generation; iGen];
                    iterateInfo.iterate       = [iterateInfo.iterate; ind];
                    iterateInfo.score         = [iterateInfo.score; genScore(ind)];
                    iterateInfo.population    = [iterateInfo.population; genPop(ind,:)];
                end
            end        
        end
        allIterateInfo{iScore} = iterateInfo;
    end
    

    if ~suppressPlot
        %Make figures
        figName = sprintf('File: %s GA Results',justFilenames{iFile});
        figure('Name',figName,'color','w')
        plot(1:length(allScore),allScore,'x','linewidth',1.5)
        grid on
        xlabel('Iterates')
        ylabel('Score')

        %typical ga fig
        figure('Name',figName,'Color','w')
        for iGen = 1:length(scoreGen)
            meanVals(iGen) = meanAfterRemovingNan(scoreGen{iGen});
            stDevVals(iGen) = std(scoreGen{iGen});
            minVals(iGen) = min(scoreGen{iGen});
            hPoint(iGen) = plot(iGen*ones(length(scoreGen{iGen}),1),scoreGen{iGen},'bx');
            hold all
        end
        h1 = plot(1:length(scoreGen),meanVals,'-bo','linewidth',1.5,'markersize',12,'markerfacecolor',0.3*ones(3,1));
        h2 = plot(1:length(scoreGen),meanVals+stDevVals,'--ro','linewidth',1.0,'markersize',12,'markerfacecolor',0.3*ones(3,1));
        lowerPoints = meanVals - stDevVals;
        for i = 1:length(scoreGen)
            if lowerPoints(i) <= minVals(i)
                lowerPoints(i) = minVals(i);
            end
        end
        
        h4 = plot(1:length(scoreGen),lowerPoints,'--ro','linewidth',1.0,'markersize',12,'markerfacecolor',0.3*ones(3,1));
        h3 = plot(1:length(scoreGen),minVals,'-go','linewidth',1.5,'markersize',12,'markerfacecolor',0.3*ones(3,1));
        
        
        clickableLegend([h1 h2 h3],'Generation Mean','Std. Dev.','Generation Best','Time Opt.')
        grid on
        xlabel('Generation')
        ylabel('Cost Function - Maneuvering time [s]')
        publicationStyleAxes
        xlim([0.9 length(scoreGen)+0.1])    
        
        axes('Position',[0.32 0.34 0.56 0.29])
        plot(1:length(scoreGen),minVals,'-go','linewidth',1.5,'markersize',12,'markerfacecolor',0.3*ones(3,1));
        hold on
        grid on
        xlabel('Generation')
        legend('Best Iterate')
        publicationStyleAxes
    %     ylabel('Cost Function - Maneuvering time [s]')
    end
end
    