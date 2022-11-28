function list = vec2list(vector,timevec)


% list = vec2list(vector,timevec)

% converts 0-1 vector of valid times (0: not valid, 1: valid)
%   into a list format of periods (i.e. [starttimes , endtimes]
% times is the vector of corresponding clock times for vector

%standardize orientation
vector = vector(:);
timevec = timevec(:);

paddedvector = [ 0 ; vector ; 0];

starttimes = timevec(find(diff(paddedvector) == 1));
endtimes = timevec(find(diff(paddedvector)==-1) - 1);

list = [starttimes , endtimes];

end