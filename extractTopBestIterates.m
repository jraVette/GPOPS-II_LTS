function extractTopBestIterates(varargin)

defaults = {'nTopIteratesToCopy',10
            'folderToCopyBestIteratesTo','bestIterates'
            'finishedRunningDirectory','finishedRunning'};
setDefaultsForVarargin(defaults,varargin);        
[~,~,~,~,~,~,~,allIterateInfo] = processGaInformationForMatlabGaV2('suppressPlot',true,'filename','gaInformation.mat');


systemCommand = sprintf('rm -r %s',folderToCopyBestIteratesTo);
system(systemCommand);
mkdir(folderToCopyBestIteratesTo);
for iIter = 1:nTopIteratesToCopy
    %Construct foldername
    for iScore = 1:length(allIterateInfo{iIter}.iterate)
        iter = allIterateInfo{iIter}.iterate(iScore);
        gen  = allIterateInfo{iIter}.generation(iScore);
        folderName = sprintf('GA_gen_%03i_iter_%03i',gen,iter);
        scoreString = strrep(num2str(allIterateInfo{iIter}.score(iScore)),'.',',');
        newName = sprintf('RANK_%03i_%s_Score_%s',iIter,folderName,scoreString);
        systemCommand = sprintf('cp -r ./%s/%s ./%s/%s',finishedRunningDirectory,folderName,folderToCopyBestIteratesTo,newName);
        system(systemCommand)
    end
    
end


