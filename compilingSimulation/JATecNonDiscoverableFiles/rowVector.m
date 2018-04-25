function row = rowVector(vector)
%This fuction will ensure the vector passed in is actually a row vector 
%
%Creation: 2014 July 18 - Jeff Anderson
%Update:   2015 July 29 - Jeff Anderson -comented out rows to allow scalers for ibox
%                         problem

[r,c] = size(vector);

% if r == 1 && c ==1 
%     error('A sclar was passed must be a vector')
if r > 1 && c > 1
    error('A matrix was passed in and must be a vector');
else
    if c > r
        row = vector';
    else
        row = vector;
    end
end