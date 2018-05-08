function extractTopBestIterates(varargin)
%Updated to use for overall or for generation
defaults = {'lookAtGeneration', []                         %Leave empty for overall, or pass in the integer value of the generation
            'nTopIteratesToCopy',10
            'folderToCopyBestIteratesTo','bestIterates'
            'finishedRunningDirectory','finishedRunning'
            'removeOldBestIterateFolder',false,
            'gaFilename','gaInformation.mat'};
setDefaultsForVarargin(defaults,varargin);        
[~,~,~,~,~,~,~,allIterateInfo] = processGaInformationForMatlabGaV2(...
    'suppressPlot',true,...
    'filename',gaFilename,...
    'processGeneration',lookAtGeneration);
load(gaFilename)
%If we want to remove all the old best iterates and check if we need to
%start a folder
if removeOldBestIterateFolder
    systemCommand = sprintf('rm -r %s',folderToCopyBestIteratesTo);
    system(systemCommand);
end

if ~exist(folderToCopyBestIteratesTo,'dir');
    mkdir(folderToCopyBestIteratesTo);
end

%Double check that nTopIteratesToCopy is less than or equal to the number
%of iterates avaiable 
if nTopIteratesToCopy > gaInfo.daq.header.populationSize; 
    nTopIteratesToCopy = gaInfo.daq.header.populationSize;
end

%Loop through iterates and copy
for iIter = 1:nTopIteratesToCopy
    %Construct foldername
    for iScore = 1:length(allIterateInfo{iIter}.iterate)
        iter = allIterateInfo{iIter}.iterate(iScore);
        gen  = allIterateInfo{iIter}.generation(iScore);
        folderName = sprintf('GA_gen_%03i_iter_%03i',gen,iter);
        scoreString = strrep(num2str(allIterateInfo{iIter}.score(iScore)),'.',',');
        if isempty(lookAtGeneration) 
            newName = sprintf('OVERALL_RANK_%03i_%s_Score_%s',iIter,folderName,scoreString);
        else
            newName = sprintf('GA_gen_%03i_RANK_%03i_iter_%03i_Score_%s',gen,iIter,iter,scoreString);
        end
        
        systemCommand = sprintf('cp -r ./%s/%s ./%s/%s',finishedRunningDirectory,folderName,folderToCopyBestIteratesTo,newName);
        system(systemCommand);
    end
    
end


