function [commonHeaders, listOfCommonChannels] = parseChannelNamesIntoMat(varargin)
%This function will parse a generic list of excel parameters in to a mat
%structure name "parameters".  This generic strucutre is fed from excel and
%very limited information needs to be known about the excel file exept that
%it must conform ot the following to parse correctly:
%
%INPUTS:
%    varargin - set defaults
%
%  - Parameters must be on the "DATA" tab
%
%  - Col A is just for formatting the template to break up
%
%  - Col B starts the paramters that will be entered in the structure
%      - If a cell in col B is empty, the row is ignored 
%      - Col B is a category label to break the parameter up into subfields 
%        in the structure to better organize them.  And the structure will
%        have a subfield to hold each item in the category
%
%  - Col C will be the paramter name - Every other field in the first row
%    will be fields in the structure and should coform to MATLAB naming
%    convention.
%
%INPUTS: NONE
%OUTPUS: 
%    commonHeaders - list of all standard channels arranged in cagegories
%    listOfCommonChannels - list of all the common channels
%
%The Excel document proivdes a template which the list to be formated like.
%
% Creation 08/02/13 - Jeff Anderson
% Updated  03/27/15 - Ryan Pawlowski - modified to parse channel names
% Updated  01/12/16 - Ryan Pawlowski - modified channelNames.xls search
%                     (line 34) to prevent common pref file errors
% Updated  04/12/16 - Jeff Anderson - had it all return a full list, non
%                     categorized
% Updated  05/18/16 - Jeff Anderson - have it saving the parsed excel data
%                     and comparing the save date to the excel file to be
%                     sure it reparses if things change in the excel.  
%                     Drastically improves performance.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

defaults  = {'userRequestReloadData',false}; %To force a reload
setDefaultsForVarargin(defaults,varargin)

if ~ispc && ~ismac || ~usejava('jvm')  %Assume this is Palmetto or headless nojvm start and just make it load the .mat file
    fprintf('Cannot parse standard xls on linux, reading .mat file of standard channels directly\n')
    load('channelNames.mat')
    reloadData = false;
else %Assumee this is normal computer
    [channelDirectory,~,~] = fileparts( mfilename('fullpath')); %Look in this directory
    fullFilename = fullfile(channelDirectory,'channelNames.xlsx');
    FileInfo = dir(fullFilename);
    xlsFileTimeStamp = FileInfo.date;
    reloadData = false;

    %Load the mat version of the file
    fullFilenameMat = fullfile(channelDirectory,'channelNames.mat');
    try 
        load(fullFilenameMat)
        if ~isequal(xlsFileTimeStamp,timeStamp)
            reloadData = true;
        end
    catch err
        reloadData = true;
    end

    %Check if the user forces a reload
    if userRequestReloadData;
        reloadData = true;
    end

end



%Re read the excel data
if reloadData
    if ~ispc
        [~,~,raw] = xlsread(fullFilename,'','','basic');
    else
        [~,~,raw] = xlsread(fullFilename);
    end

    [nRows,nCols] = size(raw);

    %Grab all the fields to enter in the structure
    % fields = cell(length(nCols),1);
    for iCol = 2:nCols
        if ~strcmp(raw{1,iCol},' ') && isempty(find(isnan(raw{1,iCol}),1))                     %See if it is a blank cell, if not enter it. 32 is an ASCII space
            fields{iCol} = fixStructureName(raw{1,iCol});                      %#ok dyn growth...not sure how big cell will be
        end
    end

    %Now enter all the data for the fields for each cahnnel
    for iRow = 2:nRows
        cellData = raw{iRow,2};
        flagNan = ~isempty(find(isnan(cellData), 1));
        flagEmpty = strcmp(cellData,' ');
        if ~flagNan  && ~flagEmpty                                                  %See if the row has parameters on it
            %Get the category
            cellData = raw{iRow,2};

            category = fixStructureName(cellData);
            category(1) = lower(category(1));

            %Get the parameter name
            cellData = raw{iRow,3};
            if isempty(find(isnan(cellData), 1)) && ~isempty(cellData) && ~strcmp(cellData,' ')
                parameter = fixStructureName(cellData);        

                %Now go through the fields
                 for iCol = 4:nCols 
                     cellData = raw{iRow,iCol};
                     if isempty(find(isnan(cellData), 1)) && ~isempty(cellData) && ~strcmp(cellData,' ')  %Make sure it's not an emptpy cell
                         if iCol >= length(fields)                                 %These are all aliases
                            commonHeaders.(category).(parameter).(fields{end}){iCol - length(fields) + 1} = cellData;
                         else                 
                            commonHeaders.(category).(parameter).(fields{iCol}) = cellData;
                         end
                     end
                 end
            end
        end    
    end
    
    %Make this into just one big list
    categories = fieldnames(commonHeaders);
    for iCat = 1:length(categories)
        channels = fieldnames(commonHeaders.(categories{iCat}));
        for iCh = 1:length(channels)
            listOfCommonChannels.(channels{iCh}) = commonHeaders.(categories{iCat}).(channels{iCh});
            listOfCommonChannels.(channels{iCh}).category = categories{iCat};
        end
    end
    
    %Save the mat file of this so we dont' have to reload each time, save
    %the exel time stamp to we can compare later
    timeStamp = xlsFileTimeStamp;
    save(fullFilenameMat,'commonHeaders','listOfCommonChannels','timeStamp')
end %Done reloading data

