function [DAQtoReturn, useGlobal] = useLocalOrGlobalDaqFiles(daq,varargin)
%This function just helps eliminate some extra code in a LOT of programs.
%Most of the daq commandline tools you can either pass the local instance
%of the daq file or let it use the globally loaded daq files. This one will
%pass along the correct one.  BE CAREFULE, if you use this, then you won't
%be able to modify the global sturcture easily!
%
%INPUTS:
%    daq - local instance of the daq file either struct or cell array
%OUTPUTS:
%    DAQtoReeturn - cell arrary of either the global or the local files
%    passed in
%    useGlobal - flag true/false weather the daqs loaded were global or
%    local
%
%Creation: 1 Oct 2015 - Jeff Anderson
%Updated:  8 Aug 2017 - Jeff Anderson - updated to include default option
%    to not loadFiles('gui') if DAQ is empty.
%Updated: 18 Aug 2017 - Jeff Anderson - option to return the selected
%    comparison DAQ if global SELECTED_DAQ obj is available
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%The way this function works, I only sometimes call global DAQ therefore,
%I'm going to turn off the warning about declarying globals before use.  I
%know what I'm doing....famous last words....
warning('off','MATLAB:declareGlobalBeforeUse')

defaults = {'promptToLoadDaq' true                                         %If non are loaded, this will prompt the user to load files using loadFiles('gui')
            'returnComparedDaqFiles',false};                               %Return the SELECTED_DAQ.CURRENT_DAQ rather than the full DAQ if global                                 
setDefaultsForVarargin(defaults,varargin);
% returnComparedDaqFiles = false;
%First see if the daq argument is ~empty, this means, we want to deal with
%a local copy:
if ~isempty(daq)
    if isa(daq,'struct')
        DAQtoReturn{1} = daq;
    elseif isa(daq,'cell')
        DAQtoReturn = daq;
    end
    useGlobal = false;
    return
end

%Now let's see if the user wanted the global comparisons
if returnComparedDaqFiles
    global SELECTED_DAQ
    if returnComparedDaqFiles && isempty(daq) && ~isempty(SELECTED_DAQ)
        daq = SELECTED_DAQ.CURRENT_DAQ;
        DAQtoReturn = daq;
        useGlobal = false;
        return
    end 
end


%First see if the global DAQ is not empty and the user wants that DAQ
global DAQ
if isempty(daq) 
    useGlobal = true; %used outside of this program
    if isempty(DAQ) && promptToLoadDaq
        status = loadFiles('gui');
        if status == 0
            error('No daq file was loaded, please load files with loadFiles()');
        end 
    end
        
    DAQtoReturn = DAQ;
    useGlobal = true;
    evalin('caller','global DAQ');
end

warning('on','MATLAB:declareGlobalBeforeUse')




    
