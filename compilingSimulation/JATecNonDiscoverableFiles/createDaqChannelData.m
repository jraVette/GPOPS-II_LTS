function channelStructure = createDaqChannelData(meas,units,name,varargin)
%This function will put the info for a typical daq entry properly in the
%daq structure.
%
%INPUTS:
%    meas - actual measurment data - should be a vector    class double nx1
%    units - string of the units for the particular channel      class char
%    name - description of channel (used in plotting)            class char
%    varargin - used to add additional fields and values see example     
%OUTPUTS:
%    channelStructure - structure to add to a daq file
%
%EXAMPLES:
%    daq.rawData.vx = createDaqChannelData([0.1 0.3 0.1],'m/s,'Velocity','mathChannel',true)
%    Output:
%        daq.rawData.vx.meas = [0.1 0.3 0.1]
%        daq.rawData.vx.units = 'm/s'
%        daq.rawData.vx.name = 'Velocity'
%        daq.rawData.vx.mathChannel = true;
%        
%Creation: 21 Oct 2014 - Jeff Anderson
%Update:   24 Nov 2014 - Jeff Anderson - to allow additional fields with
%                                        varargin

%Measurements should be vectors
[nRows,nCols] = size(meas);
if nRows > 1 && nCols > 1
    warning('daq structure designed to take vecotrs not matricies, size(meas) = %ix%i',nRows,nCols)
elseif nRows == 1 && nCols == 1
%     meas = meas;
else
    meas = rowVector(meas);
end

%Put inputs in the structure
channelStructure.meas = meas;
channelStructure.units = units;
channelStructure.name = name;

%Deal with varargin
count = 1;
if ~isempty(varargin)
    while count <= length(varargin)-1
        channelStructure.(varargin{count}) = varargin{count+1};
        count = count+2;
    end
end