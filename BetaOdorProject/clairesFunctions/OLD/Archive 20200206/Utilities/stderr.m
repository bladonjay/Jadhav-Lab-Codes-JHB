% [stderr]  = stderr(x);
%	      Returns std(x) / sqrt(length(x)-1);
%	      if x is a matrix, stderr is the standard error of each column
function [stde] = stderr(x)

stdev = std(x);
if (length(stdev) > 1)
    stde = stdev / sqrt(size(x,1)-1);
else
    stde = stdev / sqrt(length(x)-1);
end
