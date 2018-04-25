function trueFalseFlag = isStandardChannel(channelName)
%This function will see if the channel name passed in is a standard name.
%If so, return true, else false.
%INPUTS:
%    channelName - channel name                                  class char
%OUTPUS:
%    trueFalseFlag - flag either true or false if it is standard channel
%         name.                                                  class bool
%Creation: 04 Feb 2017 - Jeff Anderson
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

trueFalseFlag = false;

[~,commonChannelStruct]    = parseChannelNamesIntoMat();
commonChannels = fieldnames(commonChannelStruct);

if ismember(channelName,commonChannels)
    trueFalseFlag = true;
end