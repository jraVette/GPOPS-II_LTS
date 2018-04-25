function ind = findNearestPoint(dataVector,point)
%This funciton will find the index of the dataVector closest to the point
%of interest.
%INPUTS:
%    dataVector - input row vector of data to find point with respect to
%                                                          class double nx1  
%    point - point of interest to find the dtat point wrt  class double 1x1
%OUTPUTS
%    ind - index of this point of interest                    class int 1x1
%
%Creation: 11 Sep 2015 - Jeff Anderson
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

ind = find(min(abs(dataVector-point)) == abs(dataVector-point),1);
