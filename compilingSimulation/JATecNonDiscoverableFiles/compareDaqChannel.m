function diffRawData = compareDaqChannel(channelsToCompare,refDaq,indepVar,daqToCompare,varargin)
%This function will compare specific channels on two daq files.
%INPUTS:
%   channelsToCompare - list of channels to run comparison on 
%                                                        class char or cell
%   refDaq - reference daq instance to compare to              class struct
%   indepVar - either time or distance to serve as an interpolation basis
%                                                                class char
%   daqToCompare - (optional) local instance of a daq file to compare use
%       [] for gobal.     class struct or cell of structs (for mult. files)
%   varargin - (optional) used to set default behaviors, pass in key values
%     pairs that correspond to the defaults variable below.
%OUTPUTS:
%   diffRawData - structure that is the equivalent to the daq rawData
%       struct that constains all the differeneces calculated.
%                         class struct or cell of structs (for mult. files)
%Creation: 09 May 2017 - Jeff Anderson
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%Set defaults
defaults = {'suppressPlot',false
            'atIndepVariableTime', []};  %What range of data to look in;
setDefaultsForVarargin(defaults,varargin)

%comparison daq
if ~exist('daqToCompare','var'); daqToCompare = []; end
daqToCompare = useLocalOrGlobalDaqFiles(daqToCompare);

%Deal w/ channel inputs, if char, put in a cell
if isa(channelsToCompare,'char')
    temp = channelsToCompare;
    clear channelsToCompare;
    channelsToCompare{1} = temp;
end

%Get ref data
vars = {'tRef', indepVar};
vars = [vars; rowVector(channelsToCompare) rowVector(channelsToCompare)];
[~,refRawData] = getChannelDataFromDaqFile(refDaq,vars,'indepVariableChannel',indepVar,'atIndepVariableTime',atIndepVariableTime);



    

%Loop through files and channels and compare
for iFile = 1:length(daqToCompare)
    daq = daqToCompare{iFile};
    vars = {'indepDataChannel',indepVar};
    vars = [vars; rowVector(channelsToCompare) rowVector(channelsToCompare)];
    
    [~,rawData] = getChannelDataFromDaqFile(daq,vars,'indepVariableChannel',indepVar,'atIndepVariableTime',atIndepVariableTime);
    compRawData = rawData; %used for plotting later
    for iCh = 1:length(channelsToCompare)
        channel = channelsToCompare{iCh};
        refData = refRawData.(channel).meas;
        compData = rawData.(channel).meas;
        
        %interp to ref data
        compData = interp1(indepDataChannel,compData,tRef);
        
        %Difference
        diffCompData = refData - compData;
        
        %Put it back in the rawData struct, that'll be easy to update the
        %diffRawData struct from
        rawData.(channel).meas = diffCompData;
        rawData.(channel).rmsError = rms(diffCompData);
        rawData.(channel).integrated2NormError = trapz(tRef,diffCompData.^2);

    end
    rawData.indepDataChannel.meas = tRef;

    
    diffRawData{iFile} = rawData;
end

%If just one file, repackage
if length(diffRawData) == 1
    temp = diffRawData{1};
    clear diffRawData
    diffRawData = temp;
end


if ~suppressPlot
    for iCh = 1:length(channelsToCompare)
        channel = channelsToCompare{iCh};
        figure
        hAx = subplot(2,1,1);
        comparisonPlot.rawData = compRawData;
        plotDaqChannels('indepDataChannel',channel,comparisonPlot,'plotAxis',hAx);
        hold all
        plotRefDaq.rawData = refRawData;
        h = plotDaqChannels('tRef',channel,plotRefDaq,'plotAxis',hAx);
        set(h,'color','k')
        
        hAx = subplot(2,1,2);
        compDaq.rawData = diffRawData;
        plotDaqChannels('indepDataChannel',channel,compDaq,'plotAxis',hAx);
    end
end




