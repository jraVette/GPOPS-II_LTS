load('gaInformation.mat')
for iGen = 1:length(gaInfo.generation)
    extractTopBestIterates('lookAtGeneration',iGen)
end
extractTopBestIterates