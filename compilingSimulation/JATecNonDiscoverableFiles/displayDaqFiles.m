function varargout = displayDaqFiles(daq,varargin)
%This funciton will display the file numers and loaded daq files
%INPUTS:
%    daq (optional) - in case you want the local instance of the file
%    varargin - to set defaults
%Creation 2015 Aug 10 - Jeff Anderson
%Updated  2015 Nov 17 - Jeff Anderson - to have local file usage and
%                       suppress output to use in funciton
%Updated  2016 Aug 3 - Jeff Anderson - deal with daq and no filenames
%Updated  2017 Aug 8 - Jeff Anderson - added option to output full
%    filenames and suppress prompting to load daq if none exist.
%Updated  2017 Aug 18 - Jeff Anderson - updated options to return the
%    compared daq files as default for JATec v3.
%Updated  2018 Feb 10 - Option so it doesn't remove underscores out of the
%     label added to defaults.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%



defaults = {'suppressOutput',false
            'stripExtension',true
            'lookForShortFilename',true
            'forceCellArrayOutput',false
            'fullFileOutput',false
            'promptToLoadDaq',false %causes some issues if true, big circular loop
            'returnComparedDaqFiles',true
            'fixTextForLabel',true %Will remove underscores for plotting
            'suppressWarnings',false};
        
setDefaultsForVarargin(defaults,varargin);





if ~exist('daq','var'); 
    daq = []; 
end

needOutput = false;
if nargout>0
    needOutput = true;
end


%If user request full file output:
if fullFileOutput
    needOutput = true;
end

DAQ = useLocalOrGlobalDaqFiles(daq,'promptToLoadDaq',promptToLoadDaq,'returnComparedDaqFiles',returnComparedDaqFiles);



if ~suppressOutput
    disp(' ')
%     disp('Files:')
    disp('___________________________________________________________________')
end


filenames = cell(size(DAQ));
fullFilenames = cell(size(DAQ));

for iFile = 1:length(DAQ)
    daq = DAQ{iFile};
    
    if isfield(daq,'header')
        %Geet the filename
        if isfield(daq.header,'filename')
            filename = daq.header.filename;
        else
            if ~suppressWarnings; warning('Field ''filename'' missing, set to ''missing'''); end
            filename = 'missing';
        end
        
        %Grab the path incase we want the fullfile output
        if isfield(daq.header,'path')
            daqPath = daq.header.path;
        else
            daqPath = [];
        end
        fullFilename = fullfile(daqPath,filename);
        
        %If we want a short filename
        if lookForShortFilename
            if isfield(daq.header,'shortFilename');
                filename = daq.header.shortFilename;
            end
        end
            
    else
        if ~suppressWarnings; warning('Field ''header'' missing, cannot interpret filename'); end
        filename = 'missing';
        fullFilename = 'missing';
    end
        
    if stripExtension
        filename = strrep(filename,'.mat','');
    end

    
    if ~suppressOutput
        fprintf('File %2i: %s\n',iFile,filename);
    end
    
    if fixTextForLabel
       filename = fixLabelText(filename); 
    end
    filenames{iFile} = filename;
    fullFilenames{iFile} = fullFilename;
end

%make it a char if only one file out.
if length(filenames) == 1 && ~forceCellArrayOutput
    filenames = filenames{1};
end

if needOutput
    if fullFileOutput
        varargout{1} = fullFilenames;
    else
        varargout{1} = filenames;
    end    
end


