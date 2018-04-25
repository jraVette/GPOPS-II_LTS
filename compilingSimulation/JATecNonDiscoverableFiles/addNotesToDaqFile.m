function structOfChannelToAddNotes = addNotesToDaqFile(structOfChannelToAddNotes,notes,varargin)
%This function will add notes to the daq channel and if none are there, it
%will start them.
%INPTUS:
%    structOfChannelToAddNotes - channel structure i.e.,        class struct
%                            daq.rawData.(channels{iCh})  or,
%                            daq.header
%                        it will look for the field 'notes' in those
%                        locations
%    notes      - notes to add                                   class char
%    varargin   - used to modify defaults variable see below
%OUTPUS:
%    daqChannel - updated daq channel
%
%Creation: 22 Sep 2015 - Jeff Anderson
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

defaults = {'notesFieldName','notes'}; %This is the notes field to look for
setDefaultsForVarargin(defaults,varargin);

%If no notes exist, add that field, otherwise grab the current notes
if ~isfield(structOfChannelToAddNotes,notesFieldName)
    currentNotes = [];
else
    currentNotes = structOfChannelToAddNotes.(notesFieldName);
end

if ~isempty(currentNotes)
    notes = ['\n ' notes]; %Add the notes to a new line
end

%Concatonate
currentNotes = [currentNotes notes];


%save struct
structOfChannelToAddNotes.(notesFieldName) = currentNotes;


