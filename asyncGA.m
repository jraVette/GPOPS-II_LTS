function asyncGA(varargin)
%This function was created to use the MATLAB GA to setup new generations
%and run asyncronously.  It'll setup a generation.  Palmetto will
%evaluation. It will use that generation to setup a new one and so on until
%convergence. 
%Created: 16 Oct 2016 - Jeff Anderson (based off of previous home brewed GA
%    in jGA_remoteAsync.m.
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
addpath(genpath('./compilingSimulation'))
defaults = {'daq'       , []                %Override the daq file here
            'gaFilename','gaInformation.mat'
            'startingPointFile','startingPoint.mat'};
setDefaultsForVarargin(defaults,varargin)

genNewPopulationFlag = false;
compileScoresFlag = false;

%If we don't have this puppy started, start it up
if ~exist(gaFilename,'file')
    if isempty(daq);     
        daq = generateInitialDaq;                                          %Saves a setup of GA_information.mat
    end
    
    daq.header.gaFilename = gaFilename;                                    %This is the best spot to add that
    gaInfo.daq = daq;                                                      %Put it in teh gaInfo structure
    
    refSolution = load(sprintf('./compilingSimulation/%s',daq.header.refDaqFile));
    refSolution = refSolution.daq;
    gaInfo.referenceMeasurement = refSolution;
    
    if exist(startingPointFile,'file')                                     %See if we have a startingPoint.mat file if so, load it up
        load(startingPointFile)
        gaInfo.startingPoint = startingPoint;
    end
    
    save(gaFilename,'gaInfo')                                              %Save the file
    
end

%load the GA_information
load(gaFilename)
fprintf('Loading the GA setup from %s\n',gaFilename);
if isfield(gaInfo,'generation')
    iGen = length(gaInfo.generation);
    if ~isfield(gaInfo.generation(iGen),'Score')
        compileScoresFlag = true;
    end
else
    iGen = 1;
    genNewPopulationFlag = true;
end

%Collect scores if problem is done
if compileScoresFlag
    [gaInfo,finished] = evaluatePopulationScores(gaInfo);
    if finished
        genNewPopulationFlag = true;
    else
        fprintf('Simluations are not finished or did not hit wall time. Please re-launch or wait\n')
    end
else
    fprintf('Scores not compiled\n')
end

%Setup new generation
if genNewPopulationFlag 
    gaInfo = setupNewGeneration(gaInfo);
    fn = strrep(gaFilename,'.mat','');
    save(sprintf('%s_SnapshotAfterGen%03iSetup',fn,iGen),'gaInfo')
else
    fprintf('New generation not created\n')
end

%Save the current state of the gaInfo
save(gaFilename,'gaInfo')  ;

end%asyncGA


function gaInfo = setupNewGeneration(gaInfo)
%This function will setup the first generation randomly
%INPUTS:
%   gaInfo - the ga information structure. Specifically, it'll use the daq
%   field to make the first itration
%OUTPUT: 
%   gaInfo - updated gaInformation structure
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
daq = gaInfo.daq;

%Get the sizes of the decisions variables and bounds
c = daq.header.switchingDaq.rawData.switching.meas;
[nHorizons,nWeights] = size(c);
s = daq.header.switchingDaq.rawData.distance.meas;
lb = daq.header.switchingDaq.rawData.lb.meas;
ub = daq.header.switchingDaq.rawData.ub.meas;
nvars = numel(c);

%Constraints setup for weights
onePositions = ones(1,nWeights);
ACell = repmat({onePositions}, nHorizons, 1);
Aeq = blkdiag(ACell{:});
Beq = ones(nHorizons,1);

%Decision variables
cTime = reshape(c',[],nvars);
lb = reshape(lb',[],nvars);
ub = reshape(ub',[],nvars);


%Grab all available scores
count = 1;
c = [];
initialScores = [];

%See if there is a staringPoint field
if isfield(gaInfo,'startingPoint')
    starintPointPop = gaInfo.startingPoint.Population;
    scoreStartingPoint = gaInfo.startingPoint.scores;
    c = [c; starintPointPop];
    initialScores = [initialScores; scoreStartingPoint];
end

%Assemble all generation scores
if isfield(gaInfo,'generation')
    for iGen = 1:length(gaInfo.generation)
        genInfo = gaInfo.generation(iGen);
        cGen = gaInfo.generation(iGen).Population;
        initialScoreGen = gaInfo.generation(iGen).scores;
        c = [c; cGen];
        initialScores = [initialScores; initialScoreGen];
    end
end

%If there is still no initial guess, give it a time opt guess and empty
%score
if isempty(c)
    c = cTime;
end

%OLD WAY OF DOING THIS
% if isfield(gaInfo,'generation')                                            %If this isn't the first run, feed it the previous run
%     iGen = length(gaInfo.generation)+1;
%     lastGen = gaInfo.generation(iGen-1);
%     c = lastGen.Population;
%     initialScores = lastGen.scores;
% else
%     iGen = 1;
%     initialScores = [];
% end

%Parallel option changed on diff matlab versions
if strcmp(version('-release'),'2013a')
    useParallelFlag = 'never';
else
    useParallelFlag = false;
end

%UPDATED TO REMOVE EQUALITY CONSTRAINTS!!
Aeq = [];
Beq = [];


options = gaoptimset('InitialPopulation', c,...
                     'InitialScores',initialScores,...
                     'OutputFcn', @(options,state,flag)myOutputFunction(options,state,flag,daq),...
                     'UseParallel',useParallelFlag,...                     %So i can use a single threaded job, make this false
                     'PopulationSize',daq.header.populationSize,...
                     'Vectorized','off',...
                     'Display','diagnose',...
                     'Generations',1);
ga(@(x)fitness(x),nvars,[],[],Aeq,Beq,lb,ub,[],options);


%Needed for foldering
if isfield(gaInfo,'generation')                                            %If this isn't the first run, feed it the previous run
    iGen = length(gaInfo.generation)+1;
else
    iGen = 1;
end

%Ok, we have a population, reload the update gaInfo file and strip all the
%stuff that isn't imporant in the generation
load(daq.header.gaFilename)
tempState = gaInfo.tempState;
gaInfo.generation(iGen).Population = tempState.Population;                                       %Put the generation back
gaInfo = rmfield(gaInfo,'tempState');                                      %Remove the tempState field
save(daq.header.gaFilename,'gaInfo')                                       %Save the file

%Now, setup file all individual runs
currentGeneration = gaInfo.generation(iGen);
for iIter = 1:daq.header.populationSize
    filename = sprintf('GA_gen_%03i_iter_%03i',iGen,iIter);
    filename = fullfile('currentlyRunning',filename);
    mkdir(filename)
    sysCommand = sprintf('cp templateIterate/* %s',filename);
    system(sysCommand);
    currentDirectory = pwd;
    cd(filename);
    
    %Current decision variable
    x = currentGeneration.Population(iIter,:);
    
    %Update switching
    s0             = daq.header.initialDistance;
    controlHorizon = daq.header.controlHorizon;
    finishDistance = daq.header.finishDistance;
    s              = s0:controlHorizon:finishDistance;
    c              = reshape(x,[],length(s))';
    daq.header.switchingDaq.rawData.distance.meas  = s;
    daq.header.switchingDaq.rawData.switching.meas = c;
    
    %Save the daq file to run
    save('daqFile.mat','daq');
    
    %Go back to the right directory
    cd(currentDirectory);
end



end%setupInitialGeneration

function J = fitness(x)
%Faux fitness function so it will just generat the generation 
J = sum(x);
end

function [state,options,optchanged]  = myOutputFunction(options,state,flag,daq)

optchanged = false;

load(daq.header.gaFilename);
if isfield(gaInfo,'generation')
    iGen = length(gaInfo.generation);
else
    iGen = 1;
end

%This next part is confusion.  On the very first itration, we just passed
%in the time opt MPC guess as a partial population.  We want the rest of
%the population that MATLAB comes upwith.  So save the first time it hits
%this loop.
%On the subsequent generations however, we feed it the last population. So,
%the first time it hits this. stat.Generation == 0, it's the initial
%population we fed it.  So, wait unit generation 1.
if iGen == 1
    if state.Generation == 0
        gaInfo.tempState = state;                                           
        save(daq.header.gaFilename,'gaInfo');
    end
else
    if state.Generation == 1
        gaInfo.tempState = state;                                           
        save(daq.header.gaFilename,'gaInfo');
    end
end

end %myOutputFunction


function [updatedGaInfo,finished] = evaluatePopulationScores(gaInfo)
%This function will loop through all individuals and see if if they're done
%and extract information and store in the gaInformation structure.

%Load up the reference solution
refDaq = gaInfo.referenceMeasurement;
addpath(genpath('./compilingSimulation'));


%Assign non convergent cost after x number of minuts
daq = gaInfo.daq;
timeToGiveUpOnSim = daq.header.timeToGiveUpOnSim; %Min
updatedGaInfo = gaInfo; %If all the jobs aren't done, we don't want to update gaInfo
finished = false;

%Now, setup file all individual runs
iGen = length(gaInfo.generation);
daq = gaInfo.daq;
currentDirectory = pwd;
nonExistantSims = [];
iterToMove = [];
for iIter = 1:daq.header.populationSize
    filename = sprintf('GA_gen_%03i_iter_%03i',iGen,iIter);
    filename = fullfile('currentlyRunning',filename);
    if exist(filename,'dir')
        cd(filename)
        addpath(genpath('../../compilingSimulation'))
        if ~exist('stat.mat','file')
            fprintf('Stat file not in %s\n',filename);
            fprintf('Assume simulations have not started. Please re-start!\n')
            cd(currentDirectory)
            error('Stat file does not exist')            
        end
        load('stat.mat');
        if ~isfield(stat,'switchingDaq')                                   %this is only on finished files, so make empty field (otherwise, error will occur for mismatched structure)
            stat.switchingDaq = [];
        end
        
        %% COST FUNCTION OUTER LOOP
        %Load the solution
        solDaq = load(stat.filename);
        solDaq = solDaq.daq;
        
        if stat.simFinished && stat.conv
            diffData = compareDaqChannel({'vx';'ey'},refDaq,'distance',solDaq,'suppressPlot',true);
            stat.score = diffData.vx.integrated2NormError/1.7095e+04 + ...
                         diffData.ey.integrated2NormError/9.0544e+03;
        
        else
            stat.score = daq.header.nonConvergentCost;
        end
        
        stats(iIter) = stat;
        cd(currentDirectory)
        iterToMove = [iterToMove iIter];                                   %Used to move the simulation directory at the end
    else
        warning('Directory Does Not Exist')
        fprintf('Directory: %s Does not exist\n',filename)
        fprintf('Assume manually delted and non convergent solution\n');
        nonExistantSims = [nonExistantSims iIter];                         %Used to remove the sims 
    end
end

%Tell the nonExist sims that they didn't converge (assume manually deleted)
if ~isempty(nonExistantSims)
    for i=1:length(nonExistantSims)
        stats(nonExistantSims(i)).score = daq.header.nonConvergentCost;
        stats(nonExistantSims(i)).conv = 0;
    end
end

%Find all sims that aren't done
ind = find([stats(:).simFinished]' == 0);
cpuHrs = str2num(datestr([stats(:).elapsedTime]','HH'));
cpuMin = str2num(datestr([stats(:).elapsedTime]','MM'));
cpuSec = str2num(datestr([stats(:).elapsedTime]','SS'));
cpuTime = cpuHrs*60*60+cpuMin*60+cpuSec;
cpuTimeNotFinished = cpuTime(ind);
ind2 = find(cpuTimeNotFinished > timeToGiveUpOnSim);

%Loop through sims we need to give up on
for i = 1:length(ind2)
    iIter = ind(ind2(i));                                                  %We first found the ind of sims that were not finished, then ind2 of the ones over time
    filename = sprintf('GA_gen_%03i_iter_%03i',iGen,iIter);
    filename = fullfile('currentlyRunning',filename);
    cd(filename)    
                                                   
    cd(currentDirectory)
    stats(iIter).lapTime = daq.header.nonConvergentCost;
    stats(iIter).conv = false;
    stats(iIter).simFinished = true;
end

%Now, double check on jobs that aren't done and still haven't hit the
%cuttoff time
ind = find([stats(:).simFinished]' == 0);

%If we dont' find anymore, assume this generation is done!
if isempty(ind)
    %Ok, all sims are done let's check on the validity of scores. First
    %make sure there all min is with in a sigma of the mean.  Should be
    %within in this
    scores = [stats(:).score]';
    
    
    %Check for anomalies in scoring
%     standardDeviation = std(scores);
%     minReasonableScore = mean(scores)-standardDeviation;
%     minReasonableScore = 14.00;    
%     warningIterates = find(scores<minReasonableScore);
%     if ~isempty(warningIterates)
%         disp(' ')
%         disp('__________________________________________________________ ')
%         fprintf('Problem detected with specific iterates\n')
%         for i = 1:length(warningIterates)
%             iIter = warningIterates(i);
%             filename = sprintf('GA_gen_%03i_iter_%03i',iGen,iIter);
%             filename = fullfile('currentlyRunning',filename);
%             disp(' ')
%             fprintf('Lap time on iterate %03i = %6.4f\n',iIter,scores(iIter));
%             disp(' ')
%             fprintf('Check file: %s\n',filename);
%             disp(' ')
%             disp('__________________________________________________________ ')
%             disp(' ')
%             error('Anomalies found, check file and perhaps delete iterates')
%         end
%     end
    
    %Save results
    finished = true;
    gaInfo.generation(iGen).stats = stats';
    gaInfo.generation(iGen).scores = scores;
    updatedGaInfo = gaInfo;
    save(daq.header.gaFilename,'gaInfo')   
    
    %Clean up the files, remove the copies of all the running code
    if gaInfo.daq.header.cleanUpRunningDirecotry
        cd('currentlyRunning')
        !./cleanFinishedRuns.bsh
        cd(currentDirectory)
    end
    
    %Last thing is to move the simulation to the finishedRunningDirectory
    for i = 1:length(iterToMove)
        iIter = iterToMove(i);
        filename = sprintf('GA_gen_%03i_iter_%03i',iGen,iIter);
        filename = fullfile('currentlyRunning',filename);
        systemCommand = sprintf('mv %s finishedRunning/',filename);
        system(systemCommand);
    end
    
    %Just save the top iterates if that's what we chose to do in the
    %daqFile
    if ~isempty(gaInfo.daq.header.nIteratesToSavePerGeneration)
        extractTopBestIterates('lookAtGeneration',iGen,'nTopIteratesToCopy',gaInfo.daq.header.nIteratesToSavePerGeneration)
    end
    
    
end
    
end %evaluatePopulationScores



