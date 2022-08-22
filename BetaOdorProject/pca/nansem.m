
function sem = nansem(x,n)
% Calculate sem of data

if nargin<2, n=1; end
if n==1,
    sem = (nanstd(x,n))/sqrt(length(x)-1);
else
    sem = nanstd(x,0,n)/sqrt(size(x,n)-1);
end