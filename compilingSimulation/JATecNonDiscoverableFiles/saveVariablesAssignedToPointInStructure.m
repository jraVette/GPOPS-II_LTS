function varStructure = saveVariablesAssignedToPointInStructure(varargin)
%This funciton will save all variables assigned up to this point into a
%structure for embedding into a daq file or other result file.
%Inputs:
%    varargin - to set defaults
%
%Creation: 21 Sep 2015 - Jeff Anderson
%Updated:  02 May 2016 - Jeff Anderson - to clear variables in calling
%                        function because you shouldn't need after packing
%                        into struct
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

defaults = {'clearVariableAfterPackage',true  %clear after packaging
            'exclude',{}};  %what vars not to package
setDefaultsForVarargin(defaults,varargin)


vars = evalin('caller', 'whos');
for i = 1:length(vars)
    if ~ismember(vars(i).name,exclude)
        varStructure.(vars(i).name) = evalin('caller',vars(i).name);
        if clearVariableAfterPackage
            evalin('caller',['clear(''' vars(i).name ''')']);
        end
    end
end