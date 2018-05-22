function extractTopBestIterates(varargin)
%Updated to use for overall or for generation

defaults = {'lookAtGeneration', []                         %Leave empty for overall, or pass in the integer value of the generation
            'nTopIteratesToCopy',30
            'folderToCopyBestIteratesTo','bestIterates'
            'finishedRunningDirectory','finishedRunning'
            'removeOldBestIterateFolder',false};
setDefaultsForVarargin(defaults,varargin);        
[~,~,~,~,nonSortedScore,~,~,allIterateInfo] = processGaInformationForMatlabGaV2(...
    'suppressPlot',true,...
    'filename','gaInformation.mat',...
    'processGeneration',lookAtGeneration);

%If we want to remove all the old best iterates and check if we need to
%start a folder
if removeOldBestIterateFolder
    systemCommand = sprintf('rm -r %s',folderToCopyBestIteratesTo);
    system(systemCommand);
end

if ~exist(folderToCopyBestIteratesTo,'dir');
    mkdir(folderToCopyBestIteratesTo);
end

%Make sure that nTopIteratesToCopy isn't greater than the population
if nTopIteratesToCopy > length(allIterateInfo)
    nTopIteratesToCopy = length(allIterateInfo);
end
    

%Loop through iterates and copy
for iIter = 1:nTopIteratesToCopy
    %Construct foldername
    for iScore = 1:length(allIterateInfo{iIter}.iterate)
        iter = allIterateInfo{iIter}.iterate(iScore);
        gen  = allIterateInfo{iIter}.generation(iScore);
        searchString = sprintf('*gen_%03i*iter_%03i*',gen,iter);
        systemCommand = sprintf('find ./%s -name %s',finishedRunningDirectory,searchString);
        [~,directoryName] = system(systemCommand);
        directoryName = strsplit(directoryName);
        directoryName = directoryName{1};
        directoryName =  strtrim(directoryName);
        [~,lookForFolderWithName] = fileparts(directoryName);
%         lookForFolderWithName = sprintf('GA_gen_%03i_iter_%03i',gen,iter);
        scoreString = strrep(num2str(allIterateInfo{iIter}.score(iScore)),'.',',');
        if isempty(lookAtGeneration) 
            newName = sprintf('OVERALL_RANK_%03i_%s',iIter,lookForFolderWithName);
            if isempty(strfind(lookForFolderWithName,'Score'))
                newName = sprintf('%s_Score_%s',newName,scoreString);
            end
        else
            newName = sprintf('GA_gen_%03i_RANK_%03i_iter_%03i_Score_%s',gen,iIter,iter,scoreString);
        end
        systemCommand = sprintf('cp -r ''%s'' ''./%s/%s''',directoryName,folderToCopyBestIteratesTo,newName);
        system(systemCommand);
    end
    
end


