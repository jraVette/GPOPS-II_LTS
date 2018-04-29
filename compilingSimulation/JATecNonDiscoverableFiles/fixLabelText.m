function text = fixLabelText(input)
%This function will fix a label text by removing underscores and other
%stuff that gets messed up with latex
%
%Creation 18 July 2014 - Jeff Anderson

text = strrep(input,'_',' ');

%also, don't want to plot file extensions
% [~,justFile,~] = fileparts(text);
% text = justFile;